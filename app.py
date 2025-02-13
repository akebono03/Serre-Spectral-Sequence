import  json
from flask import Flask, render_template, request
import sqlite3
import pandas as pd
from collections import defaultdict
import re
import csv 
from io import StringIO

app = Flask(__name__)

# CSVファイルを読み込む
df_cohomology = pd.read_csv("cohomology.csv").dropna()
df_fibration = pd.read_csv("fibration.csv").dropna()
df_product = pd.read_csv("product.csv").dropna()
df_ideal = pd.read_csv("ideal.csv").dropna()

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
  generator TEXT,
  nil_exp INTEGER
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
  coe INTEGER,
  deg1 INTEGER,
  id1 INTEGER,
  deg2 INTEGER,
  id2 INTEGER,
  result INTEGER
)
""")
df_product.to_sql("product", conn, if_exists="replace", index=False)

# イデアルテーブル作成
cursor.execute("""
CREATE TABLE IF NOT EXISTS ideal (
  space TEXT,
  coe INTEGER,
  ideal_generator TEXT
)
""")
df_ideal.to_sql("ideal", conn, if_exists="replace", index=False)


conn.commit()
conn.close()

print("CSVデータをSQLiteにインポートしました。")

# def smart_split(selected_fibration):
#   # ダブルクォートで囲まれていないカンマで分割する正規表現
#   return re.split(r',(?=(?:[^"]*"[^"]*")*[^"]*$)', selected_fibration.strip())

# def smart_split(selected_fibration):
#   f = StringIO(selected_fibration)
#   reader = csv.reader(f, skipinitialspace=True)
#   return next(reader)

# def smart_split(selected_fibration):
#   res=[]
#   ok=True
#   tmp=[]
#   for c in selected_fibration:
#     if c==',':
#       if ok:
#         res.append(''.join(tmp))
#         tmp=[]
#       else:
#         tmp.append(',')
#     elif c=='{':
#       tmp.append('{')
#       ok=False
#     elif c=='}':
#       tmp.append('}')
#       ok=True
#     else:
#       tmp.append(c)
#   res.append(''.join(tmp))
#   return res

def smart_split(selected_fibration):
  s=list(selected_fibration)
  n=len(s)
  ok=True
  for i in range(n):
    if ok and s[i]==',':
      s[i] = '.'
    if s[i]=='{':
      ok=False
    if s[i]=='}':
      ok=True
  s=''.join(s)
  return s.split('.')


def get_cohomology_structure(space, type):
  conn = sqlite3.connect("cohomology.db", check_same_thread=False)
  conn.row_factory = sqlite3.Row
  cursor = conn.cursor()
  cursor.execute("SELECT deg, id, generator, nil_exp FROM cohomology WHERE space=?", (space,))
  results = cursor.fetchall()
  cursor.execute("SELECT deg1, id1, deg2, id2, result FROM product WHERE space=?", (space,))
  prod_res = cursor.fetchall()
  conn.close()

  max_deg = 30
  result_list = [[] for _ in range(max_deg + 1)]

  cohomology_dict = defaultdict(list)
  cohomology_dict[0] = ["1"]  # 0次元には1がある

  odd_generators=[]
  even_generators=[]
  gen_deg=defaultdict(int)

  for row in results:
    deg, id, nil_exp = int(row["deg"]), int(row["id"]), int(row["nil_exp"])
    a,b={"B":("a","b"), "E":("x","y"), "F":("u","v")}[type]

    generator=""
    if deg%2==1:
      if id==1:
        generator=f"{a}_{{{deg}}}"
      else:
        generator=f"{a}_{ {{deg}},{{id}} }"
    else:
      if id==1:
        generator=f"{b}_{{{deg}}}"
      else:
        generator=f"{b}_{ {{deg}},{{id}} }"

    cohomology_dict[deg].append(generator)
    gen_deg[generator]=deg
    
    # 偶数次元の生成元のi乗
    if deg % 2 == 0:
      even_generators.append(generator)
      i=2
      while deg*i<=max_deg and i<nil_exp:
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
def get_cohomology_tex(space, coefficient, type):
  conn = sqlite3.connect("cohomology.db", check_same_thread=False)
  conn.row_factory = sqlite3.Row
  cursor = conn.cursor()
  cursor.execute("SELECT id, deg, generator FROM cohomology WHERE space=? AND coe=?", (space,str(coefficient)))
  results = cursor.fetchall()
  cursor.execute("SELECT ideal_generator FROM ideal WHERE space=? AND coe=?", (space,str(coefficient)))
  ideal_res = cursor.fetchall()
  conn.close()

  if not results:
    return None, None

  if coefficient == 0:
    coe_tex = r"\mathbb{Q}"
  elif coefficient == 1:
    coe_tex = r"\mathbb{Z}"
  else:
    coe_tex = rf"\mathbb{{Z}}_{{{coefficient}}}"

  odd_generators = []
  even_generators = []
  a,b={"B":("a","b"), "E":("x","y"), "F":("u","v")}[type]

  for row in results:
    deg = int(row["deg"])
    id = int(row["id"])
    if deg % 2 == 1:
      if id == 1:
        odd_generators.append(f"{a}_{{{deg}}}")
      else:
        odd_generators.append(f"{a}_{ {{deg}}, {{id}} }")
    else:
      if id == 1:
        even_generators.append(f"{b}_{{{deg}}}")
      else:
        even_generators.append(f"{b}_{ {{deg}}, {{id}} }")

  terms = []
  ideal = []
  for row in ideal_res:
    ideal.append(row["ideal_generator"])
  if odd_generators:
    terms.append(r"\Lambda(" + ", ".join(odd_generators) + ")")
  if even_generators:
    if ideal==[]:
      terms.append(rf"{coe_tex}[" + ", ".join(even_generators) + "]")
    else:
      terms.append(rf"{coe_tex}[" + ", ".join(even_generators) + "] / " + "(" + ", ".join(ideal) + ")")

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
def get_fibration_cohomology(fibration,coefficient):
  conn = sqlite3.connect("cohomology.db", check_same_thread=False)
  conn.row_factory = sqlite3.Row
  cursor = conn.cursor()
  cursor.execute("SELECT F, E, B FROM fibration WHERE F=? AND E=? AND B=?", fibration)
  result = cursor.fetchone()
  conn.commit()
  conn.close()

  if result:
    f_coe, f_cohomology = get_cohomology_tex(result["F"],coefficient,"F")
    e_coe, e_cohomology = get_cohomology_tex(result["E"],coefficient,"E")
    b_coe, b_cohomology = get_cohomology_tex(result["B"],coefficient,"B")

    E2_tex = rf"{b_cohomology} \otimes {f_cohomology}"

    return {
      "E2": E2_tex,
      "F": (result["F"], f_coe, f_cohomology),
      "E": (result["E"], e_coe, e_cohomology),
      "B": (result["B"], b_coe, b_cohomology),
    }
  return None


def get_tensor_product(B_gens, F_gens):
  max_p, max_q = 20, 20  # グリッドサイズ
  result_grid = [[[] for _ in range(max_q)] for _ in range(max_p)]

  # 1 ⊗ 1 = 1 を設定
  result_grid[0][0].append("1")

  for p in range(max_p):
    for q in range(max_q):
      B_list = B_gens[p] if p < len(B_gens) else []
      F_list = F_gens[q] if q < len(F_gens) else []

      for b in B_list:
        for f in F_list:
          if b == "1" and f == "1":
            continue  # すでに 1 を設定済み
          elif b == "1":
            result_grid[p][q].append(f" {f}")  # 1 ⊗ f = f
          elif f == "1":
            result_grid[p][q].append(f"{b} ")  # b ⊗ 1 = b
          else:
            result_grid[p][q].append(f"{b} \otimes {f}")  # 一般形 b ⊗ f

  return result_grid


# ホームページ (ファイブレーション選択)
@app.route("/", methods=["GET", "POST"])
def index():
  fibrations = get_fibrations()
  selected_fibration = None
  selected_coefficient = "1"
  cohomologies = None

  r=request.form.get("r", "2")
  tensor_product_grid = [[[] for _ in range(20)] for _ in range(20)]

  if request.method == "POST":
    selected_fibration = request.form.get("fibration")
    selected_coefficient = request.form.get("coefficient")
    
    if selected_fibration:
      # F, E, B = selected_fibration.strip().split(",")
      F,E,B = smart_split(selected_fibration)
      cohomologies = get_fibration_cohomology((F, E, B), int(selected_coefficient))
      # print(F,E,B)
      F_gens=get_cohomology_structure(F,"F")
      B_gens=get_cohomology_structure(B,"B")
      # print(F_gens)
      # print(B_gens)

      # B_gens と F_gens は get_cohomology_structure(F), get_cohomology_structure(B) の出力
      tensor_product_grid = get_tensor_product(B_gens, F_gens)
      # for i in range(9):
      #   print(tensor_product_grid[i][:10])
      # print(tensor_product_grid[6][3])


  return render_template("index.html", fibrations=fibrations \
    , selected_fibration=selected_fibration \
    , selected_coefficient=selected_coefficient \
    , cohomologies=cohomologies \
    , tensor_product_grid=tensor_product_grid \
    , r=r)


if __name__ == "__main__":
  app.run(debug=True)
