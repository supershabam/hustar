use anyhow::Result;
#[path = "../database.rs"]
mod database;

fn main() -> Result<()> {
    let mut db = database::DatabaseMut::create("./test.bin", 4)?;
    db["accg"] += 25;
    let db = database::Database::open("./test.bin")?;
    println!("{}", db["accg"]);
    Ok(())
}
