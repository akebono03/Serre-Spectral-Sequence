from flask import Flask, render_template, request
import sqlite3
import pandas as pd
from collections import defaultdict

app = Flask(__name__)

# CSVファイルを読み込む
df_cohomology = pd.read_csv("cohomology.csv").dropna()
df_fibration = pd.read_csv("fibration.csv").dropna()
df_product = pd.read_csv("product.csv").dropna()

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

# プロダクトテーブル作成
cursor.execute("""
CREATE TABLE IF NOT EXISTS product (
  space TEXT,
  deg1 INTEGER,
  id1 INTEGER,
  deg2 INTEGER,
  id2 INTEGER,
  result INTEGER
)
""")
df_product.to_sql("product", conn, if_exists="replace", index=False)


conn.commit()
conn.close()

print("CSVデータをSQLiteにインポートしました。")


def get_cohomology_structure(space):
  conn = sqlite3.connect("cohomology.db", check_same_thread=False)
  conn.row_factory = sqlite3.Row
  cursor = conn.cursor()
  cursor.execute("SELECT deg, generator FROM cohomology WHERE space=?", (space,))
  results = cursor.fetchall()
  conn.close()

  max_deg = 30
  result_list = [[] for _ in range(max_deg + 1)]

  cohomology_dict = defaultdict(list)
  cohomology_dict[0] = ["1"]  # 0次元には1がある

  odd_generators=[]
  even_generators=[]
  gen_deg=defaultdict(int)

  for row in results:
    deg, generator = int(row["deg"]), row["generator"]
    cohomology_dict[deg].append(generator)
    gen_deg[generator]=deg
    
    # 偶数次元の生成元のi乗
    if deg % 2 == 0:
      even_generators.append(generator)
      i=2
      while deg*i<=max_deg:
        cohomology_dict[deg*i].append(f"{generator}^{i}")
        even_generators.append(f"{generator}^{i}")
        i+=1
    else:
      odd_generators.append(generator)
  
  # 生成元の積を追加
  odd_len=len(odd_generators)
  for i in range(1,1<<odd_len):
    tmp=[]
    new_deg=0
    for j in range(odd_len):
      if i>>j&1:
        tmp.append(odd_generators[j])
        new_deg+=gen_deg[odd_generators[j]]
    if len(tmp)>1:
      cohomology_dict[new_deg].append(' '.join(tmp))
  
  for deg, gens in cohomology_dict.items():
    result_list[deg] = gens
  
  return result_list

# コホモロジー環を取得する関数
def get_cohomology(space):
  conn = sqlite3.connect("cohomology.db", check_same_thread=False)
  conn.row_factory = sqlite3.Row
  cursor = conn.cursor()
  cursor.execute("SELECT coe, deg, generator FROM cohomology WHERE space=?", (space,))
  results = cursor.fetchall()
  cursor.execute("SELECT deg1, id1, deg2, id2, result FROM product WHERE space=?", (space,))
  prod_res = cursor.fetchall()
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

  # 生成元の積を追加（product.csv の情報を考慮）
  prod=[]
  for row in prod_res:
    deg1 = int(row["deg1"])
    id1 = int(row["id1"])
    deg2 = int(row["deg2"])
    id2 = int(row["id2"])
    result = int(row["result"])
    if result==0:
      if deg1==deg2 and id1==id2:
        prod.append(f"x_{{{deg1}}}^2")
      else:
        if id1==1 and id2==1:
          prod.append(f"x_{{{deg1}}} x_{{{deg2}}}")
        elif id1==1:
          prod.append(f"x_{{{deg1}}} x_{ {{deg2}},{{id2}} }")
        elif id2==1:
          prod.append(f"x_{ {{deg1}},{{id1}} } x_{{{deg2}}}")
        else:
          prod.append(f"x_{ {{deg1}},{{id1}} } x_{ {{deg2}},{{id2}} }")
    else:
      if deg1==deg2 and id1==id2:
        prod.append(f"x_{{{deg1}}}^2={{result}}")
      else:
        if id1==1 and id2==1:
          prod.append(f"x_{{{deg1}}} x_{{{deg2}}}={{result}}")
        elif id1==1:
          prod.append(f"x_{{{deg1}}} x_{ {{deg2}},{{id2}} }={{result}}")
        elif id2==1:
          prod.append(f"x_{ {{deg1}},{{id1}} } x_{{{deg2}}}={{result}}")
        else:
          prod.append(f"x_{ {{deg1}},{{id1}} } x_{ {{deg2}},{{id2}} }={{result}}")

  terms = []
  if odd_generators:
    terms.append(r"\Lambda(" + ", ".join(odd_generators) + ")")
  if even_generators:
    if prod==[]:
      terms.append(rf"{coe_tex}[" + ", ".join(even_generators) + "]")
    else:
      terms.append(rf"{coe_tex}[" + ", ".join(even_generators) + "] / " + "(" + ", ".join(prod) + ")")

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

r=2

# ホームページ (ファイブレーション選択)
@app.route("/", methods=["GET", "POST"])
def index():
  fibrations = get_fibrations()
  selected_fibration = None
  cohomologies = None

  r=request.form.get("r")

  if request.method == "POST":
    selected_fibration = request.form.get("fibration")
    if selected_fibration:
      F, E, B = selected_fibration.strip().split(",")
      cohomologies = get_fibration_cohomology((F, E, B))

      F_gens=get_cohomology_structure(F)
      B_gens=get_cohomology_structure(B)
      print(F_gens)
      print(B_gens)

  return render_template("index.html", fibrations=fibrations \
    , selected_fibration=selected_fibration, cohomologies=cohomologies \
    
    , r=r)


if __name__ == "__main__":
  app.run(debug=True)
