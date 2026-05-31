pub mod params;
pub mod cellstate;
pub mod energy;
pub mod events;
pub mod init;
pub mod grid;
pub mod dynamics;


#[macro_export]
macro_rules! trace {
    ($($args: expr),*) => {
        print!("TRACE: file: {}, line: {}", file!(), line!());
        $(
            print!(", {}: {}", stringify!($args), $args);
        )*
        println!(); // Adds a final newline
    }
}
