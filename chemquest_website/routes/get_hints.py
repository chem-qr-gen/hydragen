from chemquest_website import app, engine, meta

@app.route('/get_hints')
def get_hints():
    '''Get MS question hints from the database.'''
    hints_table = meta.tables["hints"]

    with engine.connect() as conn:
        hints = conn.execute(hints_table.select()).fetchall()
        hints = [dict(hint._mapping) for hint in hints]
    
    return hints