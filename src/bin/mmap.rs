use anyhow::Result;
#[path = "../database.rs"]
mod database;

fn main() -> Result<()> {
    let mut db = database::DatabaseMut::create("./ping.bin", 4)?;
    db["accg"] += 25;
    Ok(())
}
