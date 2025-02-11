from flask import Flask, render_template, request
import sqlite3
import pandas as pd

app = Flask(__name__)

# CSVファイルを読み込む
df_cohomology = pd.read_csv("cohomology.csv").dropna()
df_fibration = pd.read_csv("fibration.csv").dropna()

# SQLiteデータベースに接続
conn = sqlite3.connect("cohomology.db", check_same_thread=False)
cursor = conn.cursor()

# コホモロジー環テーブル作成
cursor.execute("""
CREATE TABLE IF NOT EXISTS cohomology (
  space TEXT,
  coe INTEGER,
  id INTEGER,
  deg INTEGER,
  generator TEXT
)
""")
df_cohomology.to_sql("cohomology", conn, if_exists="replace", index=False)

# ファイブレーションテーブル作成
cursor.execute("""
CREATE TABLE IF NOT EXISTS fibration (
  F TEXT,
  E TEXT,
  B TEXT,
  coe INTEGER,
  x INTEGER,
  y INTEGER,
  xid INTEGER,
  yid INTEGER,
  r INTEGER,
  imx INTEGER,
  imy INTEGER
)
""")
df_fibration.to_sql("fibration", conn, if_exists="replace", index=False)

conn.commit()
conn.close()

print("CSVデータをSQLiteにインポートしました。")


# コホモロジー環を取得する関数
def get_cohomology(space):
  conn = sqlite3.connect("cohomology.db", check_same_thread=False)
  conn.row_factory = sqlite3.Row
  cursor = conn.cursor()
  cursor.execute("SELECT coe, deg, generator FROM cohomology WHERE space=?", (space,))
  results = cursor.fetchall()
  conn.close()

  if not results:
    return None, None

  coe = results[0]["coe"]
  if coe == 0:
    coe_tex = r"\mathbb{Q}"
  elif coe == 1:
    coe_tex = r"\mathbb{Z}"
  else:
    coe_tex = rf"\mathbb{{Z}}_{{{coe}}}"

  odd_generators = []
  even_generators = []

  for row in results:
    deg = int(row["deg"])
    generator = row["generator"]
    if deg % 2 == 1:
      odd_generators.append(generator)
    else:
      even_generators.append(f"x_{{{deg}}}")

  terms = []
  if odd_generators:
    terms.append(r"\Lambda(" + ", ".join(odd_generators) + ")")
  if even_generators:
    terms.append(rf"{coe_tex}[" + ", ".join(even_generators) + "]")

  cohomology_tex = r" \otimes ".join(terms) if terms else "0"
  return coe_tex, cohomology_tex


# ファイブレーションのリストを取得
def get_fibrations():
  conn = sqlite3.connect("cohomology.db", check_same_thread=False)
  conn.row_factory = sqlite3.Row
  cursor = conn.cursor()
  cursor.execute("SELECT DISTINCT F, E, B FROM fibration")
  fibrations = cursor.fetchall()
  conn.close()
  return fibrations


# ファイブレーションに対応するコホモロジー環と E_2-term を取得
def get_fibration_cohomology(fibration):
  conn = sqlite3.connect("cohomology.db", check_same_thread=False)
  conn.row_factory = sqlite3.Row
  cursor = conn.cursor()
  cursor.execute("SELECT F, E, B FROM fibration WHERE F=? AND E=? AND B=?", fibration)
  result = cursor.fetchone()
  conn.close()

  if result:
    f_coe, f_cohomology = get_cohomology(result["F"])
    e_coe, e_cohomology = get_cohomology(result["E"])
    b_coe, b_cohomology = get_cohomology(result["B"])

    E2_tex = rf"{b_cohomology} \otimes {f_cohomology}"

    return {
      "E2": E2_tex,
      "F": (result["F"], f_coe, f_cohomology),
      "E": (result["E"], e_coe, e_cohomology),
      "B": (result["B"], b_coe, b_cohomology),
    }
  return None


# ホームページ (ファイブレーション選択)
@app.route("/", methods=["GET", "POST"])
def index():
  fibrations = get_fibrations()
  selected_fibration = None
  cohomologies = None

  if request.method == "POST":
    selected_fibration = request.form.get("fibration")
    if selected_fibration:
      F, E, B = selected_fibration.strip().split(",")
      cohomologies = get_fibration_cohomology((F, E, B))

  return render_template("index.html", fibrations=fibrations, selected_fibration=selected_fibration, cohomologies=cohomologies)


if __name__ == "__main__":
  app.run(debug=True)
