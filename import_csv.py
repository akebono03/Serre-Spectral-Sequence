import pandas as pd
import sqlite3

# CSVファイルを読み込む
df = pd.read_csv("cohomology.csv")

# 空行を削除
df = df.dropna()

# SQLiteデータベースに接続
conn = sqlite3.connect("cohomology.db")
cursor = conn.cursor()

# テーブル作成
cursor.execute("""
CREATE TABLE IF NOT EXISTS cohomology (
    space TEXT,
    coe INTEGER,
    id INTEGER,
    deg INTEGER,
    generator TEXT
)
""")

# データを挿入
df.to_sql("cohomology", conn, if_exists="replace", index=False)

# コミット＆クローズ
conn.commit()
conn.close()

print("CSVデータをSQLiteにインポートしました。")
