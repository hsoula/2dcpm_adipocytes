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
//!   Space       Reload JSON from disk
//!   Q / Esc     Quit
//!
//! Usage
//! -----
//!   cargo run --bin viewer -- data/sim3d/state_mcs000100.json

use std::sync::Arc;
use std::path::PathBuf;

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

use cpm3d::grid::SaveState;

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
    @location(1) color  : vec3<f32>,
    @location(2) normal : vec3<f32>,
};
struct VOut {
    @builtin(position) clip : vec4<f32>,
    @location(0)       col  : vec3<f32>,
};

@vertex fn vs_main(v: VIn) -> VOut {
    var out: VOut;
    out.clip = u.mvp * vec4<f32>(v.pos, 1.0);
    var col = v.color;
    if u.lighting != 0u {
        let n  = normalize(mat3x3<f32>(u.rot[0].xyz, u.rot[1].xyz, u.rot[2].xyz) * v.normal);
        let l1 = vec3<f32>(0.0, 0.0, 1.0);
        let l2 = normalize(vec3<f32>(0.6, 0.8, 0.4));
        let diff = 0.25 + 0.55 * max(dot(n, l1), 0.0) + 0.20 * max(dot(n, l2), 0.0);
        col = col * clamp(diff, 0.0, 1.0);
    }
    out.col = col;
    return out;
}

@fragment fn fs_main(in: VOut) -> @location(0) vec4<f32> {
    return vec4<f32>(in.col, 1.0);
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

#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
struct Vertex { pos: [f32;3], color: [f32;3], normal: [f32;3] }

const VERTEX_LAYOUT: wgpu::VertexBufferLayout<'static> = wgpu::VertexBufferLayout {
    array_stride: std::mem::size_of::<Vertex>() as u64,
    step_mode: wgpu::VertexStepMode::Vertex,
    attributes: &[
        wgpu::VertexAttribute { offset: 0,  shader_location: 0, format: wgpu::VertexFormat::Float32x3 },
        wgpu::VertexAttribute { offset: 12, shader_location: 1, format: wgpu::VertexFormat::Float32x3 },
        wgpu::VertexAttribute { offset: 24, shader_location: 2, format: wgpu::VertexFormat::Float32x3 },
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

fn emit_quad(out: &mut Vec<Vertex>, corners: &[[f32;3];4], col:[f32;3], nor:[f32;3], ox:f32, oy:f32, oz:f32) {
    for &i in &[0usize,1,2, 0,2,3] {
        out.push(Vertex { pos:[corners[i][0]+ox,corners[i][1]+oy,corners[i][2]+oz], color:col, normal:nor });
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
            if nb != s { emit_quad(&mut verts, corners, col, *nor, x as f32, y as f32, z as f32); }
        }
    }}}
    verts
}

fn build_centroid_mesh(grid: &[u32], w: usize, h: usize, d: usize,
                       cells: &[cpm3d::cellstate::CellState]) -> Vec<Vertex> {
    let n = cells.len();
    let mut sx = vec![0f64;n]; let mut sy = vec![0f64;n]; let mut sz = vec![0f64;n];
    let mut cnt = vec![0u64;n];
    for z in 0..d { for y in 0..h { for x in 0..w {
        let s = grid[z*w*h+y*w+x] as usize;
        if s==0||s>=n { continue; }
        sx[s]+=x as f64; sy[s]+=y as f64; sz[s]+=z as f64; cnt[s]+=1;
    }}}
    let mut verts = Vec::new();
    for (s, c) in cells.iter().enumerate() {
        if s==0||!c.alive||cnt[s]==0 { continue; }
        let (cx,cy,cz) = ((sx[s]/cnt[s] as f64) as f32,
                          (sy[s]/cnt[s] as f64) as f32,
                          (sz[s]/cnt[s] as f64) as f32);
        let r = (cnt[s] as f32).cbrt() * 0.35 + 0.5;
        let col = sigma_color(s as u32);
        for (corners, nor, _) in &CUBE_FACES {
            let sc: [[f32;3];4] = std::array::from_fn(|i|
                [corners[i][0]*2.*r+cx-r, corners[i][1]*2.*r+cy-r, corners[i][2]*2.*r+cz-r]);
            emit_quad(&mut verts, &sc, col, *nor, 0.,0.,0.);
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
    surface:        wgpu::Surface<'static>,
    surface_config: wgpu::SurfaceConfiguration,
    device:         wgpu::Device,
    queue:          wgpu::Queue,
    pipeline:       wgpu::RenderPipeline,
    depth_view:     wgpu::TextureView,
    uniform_buf:    wgpu::Buffer,
    bind_group:     wgpu::BindGroup,

    surf_buf:   Option<wgpu::Buffer>, surf_verts: u32,
    cent_buf:   Option<wgpu::Buffer>, cent_verts: u32,

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

        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label:  Some("pipeline"),
            layout: Some(&pll),
            vertex: wgpu::VertexState {
                module: &shader, entry_point: "vs_main", buffers: &[VERTEX_LAYOUT],
                compilation_options: wgpu::PipelineCompilationOptions::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader, entry_point: "fs_main",
                targets: &[Some(wgpu::ColorTargetState {
                    format: surface_config.format,
                    blend:  Some(wgpu::BlendState::REPLACE),
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
        });

        let center = Vec3::new(w as f32/2., h as f32/2., d as f32/2.);
        let dist   = w.max(h).max(d) as f32 * 2.0;
        let title  = format!("CPM3D — {} — MCS {}", json_path.display(), state.mcs);

        let mut viewer = Self {
            surface, surface_config, device, queue, pipeline, depth_view,
            uniform_buf, bind_group,
            surf_buf: None, surf_verts: 0,
            cent_buf: None, cent_verts: 0,
            camera: Camera::new(center, dist),
            drag_origin: None,
            win_size: (sz.width, sz.height),
            show_surface: true, show_centroids: false, lighting: true,
            json_path, title,
        };
        viewer.upload_meshes(&state.grid, w, h, d, &state.cells);
        viewer
    }

    fn upload_meshes(&mut self, grid:&[u32], w:usize, h:usize, d:usize,
                     cells:&[cpm3d::cellstate::CellState]) {
        let sv = build_surface_mesh(grid, w, h, d);
        let cv = build_centroid_mesh(grid, w, h, d, cells);
        self.surf_verts = sv.len() as u32;
        self.cent_verts = cv.len() as u32;
        self.surf_buf = upload(&self.device, &sv, "surf");
        self.cent_buf = upload(&self.device, &cv, "cent");
        println!("Meshes: {} surf tri, {} centroid tri",
            self.surf_verts/3, self.cent_verts/3);
    }

    fn reload(&mut self) {
        let state = load_json(&self.json_path);
        let (w,h,d) = (state.params.grid_w, state.params.grid_h, state.params.grid_d);
        self.camera.center = Vec3::new(w as f32/2., h as f32/2., d as f32/2.);
        self.title = format!("CPM3D — {} — MCS {}", self.json_path.display(), state.mcs);
        self.upload_meshes(&state.grid, w, h, d, &state.cells);
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
            if self.show_surface {
                if let Some(b) = &self.surf_buf { pass.set_vertex_buffer(0, b.slice(..)); pass.draw(0..self.surf_verts, 0..1); }
            }
            if self.show_centroids {
                if let Some(b) = &self.cent_buf { pass.set_vertex_buffer(0, b.slice(..)); pass.draw(0..self.cent_verts, 0..1); }
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
}

// ═══════════════════════════════════════════════════════════════════════════════
// Helpers
// ═══════════════════════════════════════════════════════════════════════════════

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
    println!("CPM3D Viewer  —  S: surface  C: centroids  L: lighting  R: reset  Space: reload  Q: quit");

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
