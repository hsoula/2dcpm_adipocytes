//! Interactive 3-D CPM viewer — wgpu + winit (pure Rust).
//!
//! Controls
//! --------
//!   Left-drag   Arcball rotation
//!   Scroll      Zoom
//!   S           Toggle surface
//!   C           Toggle centroids
//!   L           Toggle lighting / flat colours
//!   R           Reset camera
//!   P           Save screenshot (screenshot_NNNN.png)
//!   Space       Reload JSON from disk
//!   Q / Esc     Quit
//!
//! Usage
//! -----
//!   cargo run --bin viewer -- data/sim3d/state_mcs000100.json

use std::sync::Arc;
use std::path::PathBuf;
use png;

use clap::Parser;
use glam::{Mat4, Quat, Vec3};
use winit::{
    dpi::PhysicalPosition,
    event::{ElementState, Event, MouseButton, MouseScrollDelta, WindowEvent},
    event_loop::{ControlFlow, EventLoop},
    keyboard::{KeyCode, PhysicalKey},
    window::{Window, WindowBuilder},
};
use wgpu::util::DeviceExt;
use bytemuck::{Pod, Zeroable};

    use cpm3d::grid::{SaveState, cell_centroids};

// ═══════════════════════════════════════════════════════════════════════════════
// WGSL Shader
// ═══════════════════════════════════════════════════════════════════════════════

const SHADER_SRC: &str = r#"
struct Uniforms {
    mvp      : mat4x4<f32>,
    rot      : mat4x4<f32>,
    lighting : u32,
    _pad0    : u32,
    _pad1    : u32,
    _pad2    : u32,
};
@group(0) @binding(0) var<uniform> u: Uniforms;

struct VIn {
    @location(0) pos    : vec3<f32>,
    @location(1) color  : vec4<f32>,   // rgb + alpha (1.0 for surface, depth-faded for centroids)
    @location(2) normal : vec3<f32>,
};
struct VOut {
    @builtin(position) clip  : vec4<f32>,
    @location(0)       col   : vec3<f32>,
    @location(1)       alpha : f32,
};

@vertex fn vs_main(v: VIn) -> VOut {
    var out: VOut;
    out.clip = u.mvp * vec4<f32>(v.pos, 1.0);
    var col = v.color.rgb;
    if u.lighting != 0u {
        let n  = normalize(mat3x3<f32>(u.rot[0].xyz, u.rot[1].xyz, u.rot[2].xyz) * v.normal);
        let l1 = vec3<f32>(0.0, 0.0, 1.0);
        let l2 = normalize(vec3<f32>(0.6, 0.8, 0.4));
        let diff = 0.25 + 0.55 * max(dot(n, l1), 0.0) + 0.20 * max(dot(n, l2), 0.0);
        col = col * clamp(diff, 0.0, 1.0);
    }
    out.col   = col;
    out.alpha = v.color.a;
    return out;
}

@fragment fn fs_main(in: VOut) -> @location(0) vec4<f32> {
    return vec4<f32>(in.col, in.alpha);
}
"#;

// ═══════════════════════════════════════════════════════════════════════════════
// GPU types
// ═══════════════════════════════════════════════════════════════════════════════

#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
struct Uniforms {
    mvp:      [[f32; 4]; 4],
    rot:      [[f32; 4]; 4],
    lighting: u32,
    _pad:     [u32; 3],
}

/// pos(3) + color_rgba(4) + normal(3) = 10 floats = 40 bytes
#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
struct Vertex { pos: [f32;3], color: [f32;4], normal: [f32;3] }

const VERTEX_LAYOUT: wgpu::VertexBufferLayout<'static> = wgpu::VertexBufferLayout {
    array_stride: std::mem::size_of::<Vertex>() as u64,  // 40 bytes
    step_mode: wgpu::VertexStepMode::Vertex,
    attributes: &[
        wgpu::VertexAttribute { offset:  0, shader_location: 0, format: wgpu::VertexFormat::Float32x3 }, // pos
        wgpu::VertexAttribute { offset: 12, shader_location: 1, format: wgpu::VertexFormat::Float32x4 }, // color+alpha
        wgpu::VertexAttribute { offset: 28, shader_location: 2, format: wgpu::VertexFormat::Float32x3 }, // normal
    ],
};

// ═══════════════════════════════════════════════════════════════════════════════
// Colours
// ═══════════════════════════════════════════════════════════════════════════════

fn sigma_color(sigma: u32) -> [f32; 3] {
    if sigma == 0 { return [0.07, 0.07, 0.10]; }
    let h = (sigma as f32 * 137.508) % 360.0;
    hsv_to_rgb(h, 0.75, 0.90)
}
fn hsv_to_rgb(h: f32, s: f32, v: f32) -> [f32; 3] {
    let c = v * s;
    let x = c * (1.0 - ((h / 60.0) % 2.0 - 1.0).abs());
    let m = v - c;
    let (r,g,b) = if h<60.0{(c,x,0.0)}else if h<120.0{(x,c,0.0)}else if h<180.0{(0.0,c,x)}
                  else if h<240.0{(0.0,x,c)}else if h<300.0{(x,0.0,c)}else{(c,0.0,x)};
    [r+m, g+m, b+m]
}

// ═══════════════════════════════════════════════════════════════════════════════
// Mesh building
// ═══════════════════════════════════════════════════════════════════════════════

#[rustfmt::skip]
const CUBE_FACES: [([[f32;3];4], [f32;3], (i32,i32,i32)); 6] = [
    ([[1.,0.,0.],[1.,1.,0.],[1.,1.,1.],[1.,0.,1.]], [ 1., 0., 0.], ( 1, 0, 0)),
    ([[0.,1.,0.],[0.,0.,0.],[0.,0.,1.],[0.,1.,1.]], [-1., 0., 0.], (-1, 0, 0)),
    ([[0.,1.,0.],[1.,1.,0.],[1.,1.,1.],[0.,1.,1.]], [ 0., 1., 0.], ( 0, 1, 0)),
    ([[1.,0.,0.],[0.,0.,0.],[0.,0.,1.],[1.,0.,1.]], [ 0.,-1., 0.], ( 0,-1, 0)),
    ([[0.,0.,1.],[1.,0.,1.],[1.,1.,1.],[0.,1.,1.]], [ 0., 0., 1.], ( 0, 0, 1)),
    ([[1.,0.,0.],[0.,0.,0.],[0.,1.,0.],[1.,1.,0.]], [ 0., 0.,-1.], ( 0, 0,-1)),
];

fn emit_quad(out: &mut Vec<Vertex>, corners: &[[f32;3];4], col:[f32;3], alpha:f32, nor:[f32;3], ox:f32, oy:f32, oz:f32) {
    let color = [col[0], col[1], col[2], alpha];
    for &i in &[0usize,1,2, 0,2,3] {
        out.push(Vertex {
            pos:   [corners[i][0]+ox, corners[i][1]+oy, corners[i][2]+oz],
            color,
            normal: nor,
        });
    }
}

fn build_surface_mesh(grid: &[u32], w: usize, h: usize, d: usize) -> Vec<Vertex> {
    let mut verts = Vec::new();
    for z in 0..d { for y in 0..h { for x in 0..w {
        let s = grid[z*w*h+y*w+x];
        if s == 0 { continue; }
        let col = sigma_color(s);
        for (corners, nor, (dx,dy,dz)) in &CUBE_FACES {
            let (nx,ny,nz) = (x as i32+dx, y as i32+dy, z as i32+dz);
            let nb = if nx<0||nx>=w as i32||ny<0||ny>=h as i32||nz<0||nz>=d as i32 { 0 }
                     else { grid[nz as usize*w*h + ny as usize*w + nx as usize] };
            if nb != s { emit_quad(&mut verts, corners, col, 1.0, *nor, x as f32, y as f32, z as f32); }
        }
    }}}
    verts
}

/// Precomputed per-cell centroid info, filled on state load.
struct CentroidInfo {
    sigma:  u32,
    pos:    [f32; 3],
    radius: f32,
}

/// Precompute centroid geometry from a loaded state.
/// Calls `cell_centroids` from `grid` module — the single source of truth.
fn precompute_centroids(
    grid:  &[u32],
    w: usize, h: usize, d: usize,
    cells: &[cpm3d::cellstate::CellState],
) -> Vec<CentroidInfo> {
    let n = cells.len();
    let centroids = cell_centroids(grid, w, h, d, n);

    cells.iter().enumerate()
        .filter(|(s, c)| *s > 0 && c.alive && c.volume > 0)
        .map(|(s, c)| CentroidInfo {
            sigma:  s as u32,
            pos:    centroids[s],
            radius: 0.80, //(c.volume as f32).cbrt() * 0.35 + 0.5,
        })
        .collect()
}

/// Build the centroid cube mesh for the current frame.
///
/// Centroids are sorted back-to-front so that semi-transparent cubes blend
/// correctly.  Alpha is linearly faded from 1.0 (at `d_near`) to `ALPHA_FAR`
/// (at `d_far`) based on each centroid's distance to the camera in world space.
fn build_centroid_mesh(
    infos:    &[CentroidInfo],
    cam_pos:  Vec3,
    cam_dist: f32,
    grid_dim: f32,
) -> Vec<Vertex> {
    const ALPHA_FAR: f32 = 0.15;

    let d_near = (cam_dist - grid_dim * 0.5).max(1.0);
    let d_far  =  cam_dist + grid_dim * 0.5;

    // Collect (distance², index) for back-to-front sort
    let mut order: Vec<(f32, usize)> = infos.iter().enumerate().map(|(i, info)| {
        let p = Vec3::from(info.pos);
        ((p - cam_pos).length_squared(), i)
    }).collect();
    order.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));

    let mut verts = Vec::with_capacity(infos.len() * 36);
    for (dist_sq, i) in &order {
        let info  = &infos[*i];
        let dist  = dist_sq.sqrt();
        let t     = ((dist - d_near) / (d_far - d_near)).clamp(0.0, 1.0);
        let alpha = 1.0 - t * (1.0 - ALPHA_FAR);

        let col = sigma_color(info.sigma);
        let r   = info.radius;
        let [cx, cy, cz] = info.pos;

        for (corners, nor, _) in &CUBE_FACES {
            let sc: [[f32;3];4] = std::array::from_fn(|i|
                [corners[i][0]*2.*r+cx-r, corners[i][1]*2.*r+cy-r, corners[i][2]*2.*r+cz-r]);
            emit_quad(&mut verts, &sc, col, alpha, *nor, 0.,0.,0.);
        }
    }
    verts
}

// ═══════════════════════════════════════════════════════════════════════════════
// Camera (arcball)
// ═══════════════════════════════════════════════════════════════════════════════

struct Camera { rotation: Quat, distance: f32, center: Vec3 }

impl Camera {
    fn new(center: Vec3, distance: f32) -> Self { Self { rotation: Quat::IDENTITY, distance, center } }

    fn mvp(&self, aspect: f32) -> Mat4 {
        let proj  = Mat4::perspective_rh(45f32.to_radians(), aspect, 0.5, 5000.0);
        let view  = Mat4::from_translation(Vec3::new(0.,0.,-self.distance)) * Mat4::from_quat(self.rotation);
        let model = Mat4::from_translation(-self.center);
        proj * view * model
    }
    fn rot_mat4(&self) -> Mat4 { Mat4::from_quat(self.rotation) }

    /// Camera position in world space.
    fn world_pos(&self) -> Vec3 {
        self.center + self.rotation.inverse() * Vec3::new(0., 0., self.distance)
    }
}

fn screen_to_sphere(px:f32, py:f32, w:f32, h:f32) -> Vec3 {
    let (x,y) = (2.*px/w-1., -(2.*py/h-1.));
    let l = x*x+y*y;
    if l<=1.0 { Vec3::new(x,y,(1.-l).sqrt()) } else { Vec3::new(x,y,0.).normalize() }
}
fn arc_delta_quat(p0: Vec3, p1: Vec3) -> Quat {
    let axis = p0.cross(p1);
    if axis.length_squared() < 1e-12 { return Quat::IDENTITY; }
    Quat::from_axis_angle(axis.normalize(), p0.dot(p1).clamp(-1.,1.).acos())
}

// ═══════════════════════════════════════════════════════════════════════════════
// Viewer state
// ═══════════════════════════════════════════════════════════════════════════════

struct Viewer {
    surface:             wgpu::Surface<'static>,
    surface_config:      wgpu::SurfaceConfiguration,
    device:              wgpu::Device,
    queue:               wgpu::Queue,
    pipeline:            wgpu::RenderPipeline,
    screenshot_pipeline: wgpu::RenderPipeline,
    screenshot_counter:  u32,
    depth_view:          wgpu::TextureView,
    uniform_buf:         wgpu::Buffer,
    bind_group:          wgpu::BindGroup,

    surf_buf:   Option<wgpu::Buffer>, surf_verts: u32,

    /// Precomputed centroid positions/radii — rebuilt only on state reload.
    centroid_data: Vec<CentroidInfo>,
    /// Largest grid dimension — used to compute depth-fade range.
    grid_dim: f32,

    camera:       Camera,
    drag_origin:  Option<(f32,f32)>,
    win_size:     (u32,u32),

    show_surface:   bool,
    show_centroids: bool,
    lighting:       bool,

    json_path: PathBuf,
    title:     String,
}

impl Viewer {
    async fn new(window: Arc<Window>, json_path: PathBuf) -> Self {
        let state = load_json(&json_path);
        let (w,h,d) = (state.params.grid_w, state.params.grid_h, state.params.grid_d);

        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::default());
        let surface  = instance.create_surface(window.clone()).unwrap();

        let adapter = instance.request_adapter(&wgpu::RequestAdapterOptions {
            compatible_surface: Some(&surface),
            power_preference: wgpu::PowerPreference::HighPerformance,
            force_fallback_adapter: false,
        }).await.expect("no GPU adapter found");

        let (device, queue) = adapter.request_device(&wgpu::DeviceDescriptor::default(), None)
            .await.expect("device request failed");

        let sz   = window.inner_size();
        let caps = surface.get_capabilities(&adapter);
        let fmt  = caps.formats.iter().copied().find(|f| f.is_srgb()).unwrap_or(caps.formats[0]);

        let surface_config = wgpu::SurfaceConfiguration {
            usage:        wgpu::TextureUsages::RENDER_ATTACHMENT,
            format:       fmt,
            width:        sz.width,
            height:       sz.height,
            present_mode: wgpu::PresentMode::Fifo,
            alpha_mode:   caps.alpha_modes[0],
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(&device, &surface_config);

        let depth_view = make_depth_view(&device, sz.width, sz.height);

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label:  Some("shader"),
            source: wgpu::ShaderSource::Wgsl(SHADER_SRC.into()),
        });

        let uniform_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("uniforms"),
            size:  std::mem::size_of::<Uniforms>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("bgl"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0, visibility: wgpu::ShaderStages::VERTEX,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false, min_binding_size: None,
                },
                count: None,
            }],
        });
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("bg"), layout: &bgl,
            entries: &[wgpu::BindGroupEntry { binding: 0, resource: uniform_buf.as_entire_binding() }],
        });

        let pll = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("pll"), bind_group_layouts: &[&bgl], push_constant_ranges: &[],
        });

        let pipeline            = make_pipeline(&device, &shader, &pll, surface_config.format);
        let screenshot_pipeline = make_pipeline(&device, &shader, &pll, wgpu::TextureFormat::Rgba8Unorm);

        let center   = Vec3::new(w as f32/2., h as f32/2., d as f32/2.);
        let dist     = w.max(h).max(d) as f32 * 2.0;
        let grid_dim = w.max(h).max(d) as f32;
        let title    = format!("CPM3D — {} — MCS {}", json_path.display(), state.mcs);

        let centroid_data = precompute_centroids(&state.grid, w, h, d, &state.cells);
        let sv = build_surface_mesh(&state.grid, w, h, d);

        let surf_verts = sv.len() as u32;
        let surf_buf   = upload(&device, &sv, "surf");

        println!("Meshes: {} surf tri, {} centroids", surf_verts/3, centroid_data.len());

        Self {
            surface, surface_config, device, queue, pipeline, screenshot_pipeline,
            screenshot_counter: 0,
            depth_view, uniform_buf, bind_group,
            surf_buf, surf_verts,
            centroid_data, grid_dim,
            camera: Camera::new(center, dist),
            drag_origin: None,
            win_size: (sz.width, sz.height),
            show_surface: true, show_centroids: false, lighting: true,
            json_path, title,
        }
    }

    fn reload(&mut self) {
        let state = load_json(&self.json_path);
        let (w,h,d) = (state.params.grid_w, state.params.grid_h, state.params.grid_d);
        self.camera.center = Vec3::new(w as f32/2., h as f32/2., d as f32/2.);
        self.grid_dim = w.max(h).max(d) as f32;
        self.title = format!("CPM3D — {} — MCS {}", self.json_path.display(), state.mcs);

        let sv = build_surface_mesh(&state.grid, w, h, d);
        self.surf_verts = sv.len() as u32;
        self.surf_buf   = upload(&self.device, &sv, "surf");
        self.centroid_data = precompute_centroids(&state.grid, w, h, d, &state.cells);

        println!("Reloaded: {} surf tri, {} centroids", self.surf_verts/3, self.centroid_data.len());
    }

    fn resize(&mut self, nw: u32, nh: u32) {
        if nw == 0 || nh == 0 { return; }
        self.surface_config.width  = nw;
        self.surface_config.height = nh;
        self.surface.configure(&self.device, &self.surface_config);
        self.depth_view = make_depth_view(&self.device, nw, nh);
        self.win_size = (nw, nh);
    }

    fn render(&mut self, window: &Window) {
        let aspect = self.win_size.0 as f32 / self.win_size.1 as f32;
        let uni = Uniforms {
            mvp: self.camera.mvp(aspect).to_cols_array_2d(),
            rot: self.camera.rot_mat4().to_cols_array_2d(),
            lighting: self.lighting as u32,
            _pad: [0;3],
        };
        self.queue.write_buffer(&self.uniform_buf, 0, bytemuck::bytes_of(&uni));

        // Build centroid buffer for this frame (cheap — few dozen cubes).
        // Rebuilt every frame so alpha updates as camera rotates.
        let (cent_buf, cent_verts) = if self.show_centroids && !self.centroid_data.is_empty() {
            let cv = build_centroid_mesh(
                &self.centroid_data,
                self.camera.world_pos(),
                self.camera.distance,
                self.grid_dim,
            );
            let n = cv.len() as u32;
            (upload(&self.device, &cv, "cent"), n)
        } else {
            (None, 0u32)
        };

        let output = match self.surface.get_current_texture() {
            Ok(t)  => t,
            Err(_) => { self.surface.configure(&self.device, &self.surface_config); return; }
        };
        let view = output.texture.create_view(&wgpu::TextureViewDescriptor::default());
        let mut enc = self.device.create_command_encoder(&wgpu::CommandEncoderDescriptor::default());
        {
            let mut pass = enc.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: None,
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view, resolve_target: None,
                    ops: wgpu::Operations {
                        load:  wgpu::LoadOp::Clear(wgpu::Color { r:0.06, g:0.06, b:0.09, a:1.0 }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &self.depth_view,
                    depth_ops:   Some(wgpu::Operations { load: wgpu::LoadOp::Clear(1.0), store: wgpu::StoreOp::Store }),
                    stencil_ops: None,
                }),
                ..Default::default()
            });
            pass.set_pipeline(&self.pipeline);
            pass.set_bind_group(0, &self.bind_group, &[]);

            // Draw opaque surface first, then transparent centroids on top.
            if self.show_surface {
                if let Some(b) = &self.surf_buf {
                    pass.set_vertex_buffer(0, b.slice(..));
                    pass.draw(0..self.surf_verts, 0..1);
                }
            }
            if self.show_centroids {
                if let Some(b) = &cent_buf {
                    pass.set_vertex_buffer(0, b.slice(..));
                    pass.draw(0..cent_verts, 0..1);
                }
            }
        }
        self.queue.submit(std::iter::once(enc.finish()));
        output.present();

        let mode = match (self.show_surface, self.show_centroids) {
            (true,  true)  => "surface+centroids",
            (true,  false) => "surface",
            (false, true)  => "centroids",
            _              => "—",
        };
        window.set_title(&format!("{} [{mode}] [{}]", self.title, if self.lighting {"lit"} else {"flat"}));
    }

    fn mouse_press(&mut self, p: PhysicalPosition<f64>) {
        self.drag_origin = Some((p.x as f32, p.y as f32));
    }
    fn mouse_release(&mut self) { self.drag_origin = None; }
    fn mouse_move(&mut self, p: PhysicalPosition<f64>) {
        if let Some((ox,oy)) = self.drag_origin {
            let (w,h) = (self.win_size.0 as f32, self.win_size.1 as f32);
            let delta = arc_delta_quat(screen_to_sphere(ox,oy,w,h), screen_to_sphere(p.x as f32,p.y as f32,w,h));
            self.camera.rotation = (delta * self.camera.rotation).normalize();
            self.drag_origin = Some((p.x as f32, p.y as f32));
        }
    }
    fn scroll(&mut self, d: f32) {
        self.camera.distance = (self.camera.distance * (1.0 - d * 0.1)).max(1.0);
    }

    fn screenshot(&mut self) {
        let (w, h) = self.win_size;

        let tex = self.device.create_texture(&wgpu::TextureDescriptor {
            label: Some("ss_tex"),
            size:  wgpu::Extent3d { width: w, height: h, depth_or_array_layers: 1 },
            mip_level_count: 1, sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format:    wgpu::TextureFormat::Rgba8Unorm,
            usage:     wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
            view_formats: &[],
        });
        let tex_view   = tex.create_view(&wgpu::TextureViewDescriptor::default());
        let depth_view = make_depth_view(&self.device, w, h);

        let aspect = w as f32 / h as f32;
        let uni = Uniforms {
            mvp:      self.camera.mvp(aspect).to_cols_array_2d(),
            rot:      self.camera.rot_mat4().to_cols_array_2d(),
            lighting: self.lighting as u32,
            _pad:     [0; 3],
        };
        self.queue.write_buffer(&self.uniform_buf, 0, bytemuck::bytes_of(&uni));

        let (cent_buf, cent_verts) = if self.show_centroids && !self.centroid_data.is_empty() {
            let cv = build_centroid_mesh(
                &self.centroid_data, self.camera.world_pos(),
                self.camera.distance, self.grid_dim,
            );
            let n = cv.len() as u32;
            (upload(&self.device, &cv, "ss_cent"), n)
        } else {
            (None, 0u32)
        };

        let mut enc = self.device.create_command_encoder(
            &wgpu::CommandEncoderDescriptor { label: Some("ss_enc") });
        {
            let mut pass = enc.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("ss_pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &tex_view, resolve_target: None,
                    ops: wgpu::Operations {
                        load:  wgpu::LoadOp::Clear(wgpu::Color { r: 0.06, g: 0.06, b: 0.09, a: 1.0 }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0), store: wgpu::StoreOp::Store }),
                    stencil_ops: None,
                }),
                ..Default::default()
            });
            pass.set_pipeline(&self.screenshot_pipeline);
            pass.set_bind_group(0, &self.bind_group, &[]);
            if self.show_surface {
                if let Some(b) = &self.surf_buf {
                    pass.set_vertex_buffer(0, b.slice(..));
                    pass.draw(0..self.surf_verts, 0..1);
                }
            }
            if self.show_centroids {
                if let Some(b) = &cent_buf {
                    pass.set_vertex_buffer(0, b.slice(..));
                    pass.draw(0..cent_verts, 0..1);
                }
            }
        }

        // Copy to CPU-readable staging buffer (rows must be 256-byte aligned)
        let row_bytes         = w * 4;
        let align             = wgpu::COPY_BYTES_PER_ROW_ALIGNMENT;
        let row_bytes_aligned = (row_bytes + align - 1) / align * align;
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("ss_staging"),
            size:  (row_bytes_aligned * h) as u64,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        enc.copy_texture_to_buffer(
            tex.as_image_copy(),
            wgpu::ImageCopyBuffer {
                buffer: &staging,
                layout: wgpu::ImageDataLayout {
                    offset:         0,
                    bytes_per_row:  Some(row_bytes_aligned),
                    rows_per_image: Some(h),
                },
            },
            wgpu::Extent3d { width: w, height: h, depth_or_array_layers: 1 },
        );
        self.queue.submit(std::iter::once(enc.finish()));

        let slice = staging.slice(..);
        slice.map_async(wgpu::MapMode::Read, |_| {});
        self.device.poll(wgpu::Maintain::Wait);

        let raw = slice.get_mapped_range();
        let mut pixels: Vec<u8> = Vec::with_capacity((w * h * 4) as usize);
        for row in 0..h as usize {
            let start = row * row_bytes_aligned as usize;
            pixels.extend_from_slice(&raw[start..start + row_bytes as usize]);
        }
        drop(raw);
        staging.unmap();

        let path = format!("screenshot_{:04}.png", self.screenshot_counter);
        self.screenshot_counter += 1;
        let file = std::fs::File::create(&path).expect("cannot create screenshot file");
        let mut encoder = png::Encoder::new(file, w, h);
        encoder.set_color(png::ColorType::Rgba);
        encoder.set_depth(png::BitDepth::Eight);
        encoder.write_header().unwrap().write_image_data(&pixels).unwrap();
        println!("Screenshot → {path}");
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Helpers
// ═══════════════════════════════════════════════════════════════════════════════

fn make_pipeline(
    device: &wgpu::Device,
    shader: &wgpu::ShaderModule,
    pll:    &wgpu::PipelineLayout,
    fmt:    wgpu::TextureFormat,
) -> wgpu::RenderPipeline {
    device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label:  Some("pipeline"),
        layout: Some(pll),
        vertex: wgpu::VertexState {
            module: shader, entry_point: "vs_main", buffers: &[VERTEX_LAYOUT],
            compilation_options: wgpu::PipelineCompilationOptions::default(),
        },
        fragment: Some(wgpu::FragmentState {
            module: shader, entry_point: "fs_main",
            targets: &[Some(wgpu::ColorTargetState {
                format: fmt,
                blend:  Some(wgpu::BlendState::ALPHA_BLENDING),
                write_mask: wgpu::ColorWrites::ALL,
            })],
            compilation_options: wgpu::PipelineCompilationOptions::default(),
        }),
        primitive: wgpu::PrimitiveState {
            topology: wgpu::PrimitiveTopology::TriangleList,
            cull_mode: None,
            ..Default::default()
        },
        depth_stencil: Some(wgpu::DepthStencilState {
            format: wgpu::TextureFormat::Depth32Float,
            depth_write_enabled: true,
            depth_compare: wgpu::CompareFunction::Less,
            stencil: wgpu::StencilState::default(),
            bias: wgpu::DepthBiasState::default(),
        }),
        multisample: wgpu::MultisampleState::default(),
        multiview: None,
    })
}

fn make_depth_view(device: &wgpu::Device, w: u32, h: u32) -> wgpu::TextureView {
    device.create_texture(&wgpu::TextureDescriptor {
        label: Some("depth"),
        size:  wgpu::Extent3d { width: w, height: h, depth_or_array_layers: 1 },
        mip_level_count: 1, sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format:    wgpu::TextureFormat::Depth32Float,
        usage:     wgpu::TextureUsages::RENDER_ATTACHMENT,
        view_formats: &[],
    }).create_view(&wgpu::TextureViewDescriptor::default())
}

fn upload(device: &wgpu::Device, verts: &[Vertex], label: &'static str) -> Option<wgpu::Buffer> {
    if verts.is_empty() { return None; }
    Some(device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some(label), contents: bytemuck::cast_slice(verts), usage: wgpu::BufferUsages::VERTEX,
    }))
}

fn load_json(path: &PathBuf) -> SaveState {
    let json = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
    serde_json::from_str(&json)
        .unwrap_or_else(|e| panic!("cannot parse {}: {e}", path.display()))
}

// ═══════════════════════════════════════════════════════════════════════════════
// CLI + main
// ═══════════════════════════════════════════════════════════════════════════════

#[derive(Parser)]
struct Cli {
    /// Path to a cpm3d JSON snapshot
    json: PathBuf,
}

fn main() {
    let cli = Cli::parse();
    println!("CPM3D Viewer  —  S: surface  C: centroids  L: lighting  R: reset  P: screenshot  Space: reload  Q: quit");

    let event_loop = EventLoop::new().unwrap();
    let window = Arc::new(
        WindowBuilder::new()
            .with_title("CPM3D Viewer")
            .with_inner_size(winit::dpi::LogicalSize::new(1024u32, 768u32))
            .build(&event_loop)
            .unwrap(),
    );

    let mut viewer = pollster::block_on(Viewer::new(window.clone(), cli.json));
    let mut cursor  = PhysicalPosition::new(0.0f64, 0.0f64);

    event_loop.run(move |event, target| {
        target.set_control_flow(ControlFlow::Poll);

        match event {
            Event::WindowEvent { event: we, .. } => match we {
                WindowEvent::CloseRequested => target.exit(),
                WindowEvent::RedrawRequested => viewer.render(&window),
                WindowEvent::Resized(sz)     => viewer.resize(sz.width, sz.height),

                WindowEvent::CursorMoved { position, .. } => {
                    cursor = position;
                    viewer.mouse_move(position);
                }
                WindowEvent::MouseInput { state: ElementState::Pressed,  button: MouseButton::Left, .. } => viewer.mouse_press(cursor),
                WindowEvent::MouseInput { state: ElementState::Released, button: MouseButton::Left, .. } => viewer.mouse_release(),
                WindowEvent::MouseWheel { delta, .. } => {
                    let d = match delta {
                        MouseScrollDelta::LineDelta(_, y)  => y,
                        MouseScrollDelta::PixelDelta(p) => p.y as f32 * 0.01,
                    };
                    viewer.scroll(d);
                }

                WindowEvent::KeyboardInput {
                    event: winit::event::KeyEvent {
                        physical_key: PhysicalKey::Code(key),
                        state: ElementState::Pressed, ..
                    }, ..
                } => match key {
                    KeyCode::Escape | KeyCode::KeyQ => target.exit(),
                    KeyCode::KeyS => {
                        viewer.show_surface = !viewer.show_surface;
                        if !viewer.show_surface && !viewer.show_centroids { viewer.show_centroids = true; }
                        println!("Surface: {}  Centroids: {}", viewer.show_surface, viewer.show_centroids);
                    }
                    KeyCode::KeyC => {
                        viewer.show_centroids = !viewer.show_centroids;
                        if !viewer.show_surface && !viewer.show_centroids { viewer.show_surface = true; }
                        println!("Surface: {}  Centroids: {}", viewer.show_surface, viewer.show_centroids);
                    }
                    KeyCode::KeyL => { viewer.lighting = !viewer.lighting; println!("Lighting: {}", viewer.lighting); }
                    KeyCode::KeyR => { viewer.camera.rotation = Quat::IDENTITY; println!("Camera reset"); }
                    KeyCode::KeyP => { viewer.screenshot(); }
                    KeyCode::Space => { println!("Reloading…"); viewer.reload(); }
                    _ => {}
                },
                _ => {}
            },
            Event::AboutToWait => window.request_redraw(),
            _ => {}
        }
    }).unwrap();
}
