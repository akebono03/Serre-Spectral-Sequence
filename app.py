import  json
from flask import Flask, render_template, request
import sqlite3
import pandas as pd
from collections import defaultdict
import re
import csv 
from io import StringIO
from itertools import product

app = Flask(__name__)

# CSVファイルを読み込む
df_cohomology = pd.read_csv("cohomology.csv").dropna()
df_fibration = pd.read_csv("fibration.csv").dropna()
df_ideal = pd.read_csv("ideal.csv").dropna()
df_reference = pd.read_csv("reference.csv").dropna()

# SQLiteデータベースに接続
conn = sqlite3.connect("cohomology.db", check_same_thread=False)
cursor = conn.cursor()

# コホモロジー環テーブル作成
cursor.execute("""
CREATE TABLE IF NOT EXISTS cohomology (
  space TEXT,
  coe INTEGER,
  type INTEGER,
  deg INTEGER,
  id INTEGER,
  generator TEXT,
  nil_exp INTEGER,
  gen_order INTEGER
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
  r INTEGER,
  kill TEXT
)
""")
df_fibration.to_sql("fibration", conn, if_exists="replace", index=False)

# イデアルテーブル作成
cursor.execute("""
CREATE TABLE IF NOT EXISTS ideal (
  space TEXT,
  coe INTEGER,
  ideal_generator TEXT
)
""")
df_ideal.to_sql("ideal", conn, if_exists="replace", index=False)

# referenceテーブル作成
cursor.execute("""
CREATE TABLE IF NOT EXISTS reference (
  F TEXT,
  E TEXT,
  B TEXT,
  coe INTEGER,
  reference TEXT
)
""")
df_reference.to_sql("reference", conn, if_exists="replace", index=False)


conn.commit()
conn.close()

print("CSVデータをSQLiteにインポートしました。")

max_deg = 30


def smart_split(selected_fibration):
  s=list(selected_fibration)
  n=len(s)
  ok=True
  for i in range(n):
    if ok and s[i]==',':
      s[i] = '.'
    if s[i]=='{' or s[i]=='(':
      ok=False
    if s[i]=='}' or s[i]==')':
      ok=True
  s=''.join(s)
  return s.split('.')

def get_generator_list(space, coefficient, space_type):
  conn = sqlite3.connect("cohomology.db", check_same_thread=False)
  conn.row_factory = sqlite3.Row
  cursor = conn.cursor()
  cursor.execute("SELECT gen_type, deg, id, nil_exp, g_order FROM cohomology WHERE space=? AND coe=?", (space,coefficient,))
  results = cursor.fetchall()
  conn.close()

  id_list = [[] for _ in range(max_deg+1)]
  gen_list = []
  gen_deg = defaultdict(int)
  gen_order = defaultdict(int)
  gen_nil = defaultdict(int)
  is_except = False

  for row in results:
    gen_type, deg, id, nil_exp, order = int(row["gen_type"]), int(row["deg"]), int(row["id"]), int(row["nil_exp"]), int(row["g_order"])
    id_list[deg].append((id,nil_exp,order))
    if gen_type == 5:
      is_except = True

  deg_only_one = True
  cnt=0
  for deg in range(max_deg+1):
    if len(id_list[deg])>0:
      cnt+=1
  if cnt>1:
    deg_only_one = False

  for deg in range(max_deg+1):
    if len(id_list[deg])==0: continue
    x = {"B":("b","a"), "E":("y","x"), "F":("v","u")}[space_type][deg%2]
    if len(id_list[deg])==1:
      gen = x if deg_only_one else f"{x}_{{{deg}}}"
      gen_list.append(gen)
      gen_deg[gen] = deg
      _, nil_exp, order = id_list[deg][0]
      gen_order[gen] = order
      gen_nil[gen] = nil_exp
    else:
      for id, nil_exp, order in id_list[deg]:
        gen = f"{x}_{{{id}}}" if deg_only_one else x + "_{" + str(deg) + "," + str(id) + "}"
        gen_list.append(gen)
        gen_deg[gen] = deg
        gen_order[gen] = order
        gen_nil[gen] = nil_exp

  return gen_list, gen_deg, gen_order, gen_nil, is_except


def get_cohomology_structure(space, coefficient, space_type):
  conn = sqlite3.connect("cohomology.db", check_same_thread=False)
  conn.row_factory = sqlite3.Row
  cursor = conn.cursor()
  cursor.execute("SELECT gen_type, deg, id, generator, nil_exp FROM cohomology WHERE space=? AND coe=?", (space,coefficient,))
  results = cursor.fetchall()
  conn.close()

  result_list = [[] for _ in range(max_deg + 1)]

  cohomology_dict = defaultdict(list)
  cohomology_dict[0] = ["1"]  # 0次元には1がある

  odd_generators=[]
  even_generators=[]

  gen_list, gen_deg, _, gen_nil, is_except = get_generator_list(space, coefficient, space_type)

  # 例外の場合
  if is_except:
    for gen in gen_list:
      result_list[gen_deg[gen]].append(gen)
    result_list[0].append('1')
    return result_list

  # 通常の場合
  for gen in gen_list:
    if gen_deg[gen]%2 == 1:
      odd_generators.append(gen)
    else:
      even_generators.append(gen)

  nil_list = []
  for odd_gen in odd_generators:
    nil_list.append([1,odd_gen])
  for even_gen in even_generators:
    nil_list.append([1, even_gen])
    i=2
    deg = gen_deg[even_gen]
    while deg*i<=max_deg and i<gen_nil[even_gen]:
      nil_list[-1].append(f"{even_gen}^{{{i}}}")
      gen_deg[f"{even_gen}^{{{i}}}"] = deg*i
      i+=1

  combinations = list(product(*nil_list))
  for comb in combinations:
    new_deg = 0
    new_el = []
    for el in comb:
      if el == 1: continue
      new_deg += gen_deg[el]
      new_el.append(el)
    if len(new_el)>0:
      cohomology_dict[new_deg].append(" ".join(new_el))

  for deg in cohomology_dict.keys():
    gens = cohomology_dict[deg]
    gens.sort()
    if deg <= max_deg:
      result_list[deg] = gens

  result_list[0]=['1'] # 必要
  return result_list

def get_ideal(space, coefficient, space_type):
  conn = sqlite3.connect("cohomology.db", check_same_thread=False)
  conn.row_factory = sqlite3.Row
  cursor = conn.cursor()
  # cursor.execute("SELECT gen_type, deg, id, generator, g_order FROM cohomology WHERE space=? AND coe=?", (space,str(coefficient)))
  # results = cursor.fetchall()
  cursor.execute("SELECT ideal_generator FROM ideal WHERE space=? AND coe=?", (space,str(coefficient)))
  ideal_res = cursor.fetchall()
  conn.close()

  ideal = []
  for row in ideal_res:
    ideal_gen = row["ideal_generator"]
    x,y = {"B":("b","a"), "E":("y","x"), "F":("v","u")}[space_type]
    tmp = []
    for c in ideal_gen:
      if c in {"b","y","v"}:
        tmp.append(x)
      elif c in {"a","x","u"}:
        tmp.append(y)
      else:
        tmp.append(c) 
    new_ideal_gen = "".join(tmp)
    ideal.append(new_ideal_gen)
  return ideal


# コホモロジー環の係数と環構造を取得する関数
def get_cohomology_tex(space, coefficient, space_type):
  conn = sqlite3.connect("cohomology.db", check_same_thread=False)
  conn.row_factory = sqlite3.Row
  cursor = conn.cursor()
  cursor.execute("SELECT gen_type, deg, id, generator, g_order FROM cohomology WHERE space=? AND coe=?", (space,str(coefficient)))
  results = cursor.fetchall()
  conn.close()

  if not results:
    return None, None

  if coefficient == 0:
    coe_tex = r"\mathbb{Q}"
  elif coefficient == 1:
    coe_tex = r"\mathbb{Z}"
  else:
    coe_tex = rf"\mathbb{{Z}}_{{{coefficient}}}"

  # 一点の場合
  if space == '*':
    return coe_tex, coe_tex

  gen_list, gen_deg, gen_order, _, is_except = get_generator_list(space, coefficient, space_type)

  # 例外の場合
  if is_except:
    groups = []
    for gen in gen_list:
      if gen_order[gen] == 1:
        groups.append("Z\{" + gen + "\}")
      else:
        groups.append(f"Z_{{{gen_order[gen]}}}" + "\{" + gen + "\}")      

    groups.append("\cdots")
    cohomology_tex = r" \oplus ".join(groups) if groups else "0"
    return coe_tex, "Z\{1\}\\oplus " + cohomology_tex

  # 通常の場合
  odd_generators = []
  even_generators = []
  terms = []

  for gen in gen_list:
    if gen_deg[gen]%2 == 1:
      odd_generators.append(gen)
    else:
      even_generators.append(gen)

  ideal = get_ideal(space, coefficient, space_type)
  
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

def str_to_tuple(s):
  """ '(a,b)' のような文字列を (a, b) のタプルに変換する """
  s = s.strip("()")  # 括弧を除去
  elements = s.split(",")  # カンマで分割
  parsed_elements = [int(e) if e.strip().isdigit() else e.strip() for e in elements]  # 型変換
  return tuple(parsed_elements)

def get_deleted(F,E,B,coe):
  conn = sqlite3.connect("cohomology.db", check_same_thread=False)
  conn.row_factory = sqlite3.Row
  cursor = conn.cursor()
  cursor.execute("SELECT r,kill FROM fibration WHERE F=? AND E=? AND B=? AND coe=?", (F,E,B, coe))
  result = cursor.fetchall()
  conn.commit()
  conn.close()

  dels = [set() for _ in range(max_deg+1)]

  for row in result:
    dels[row[0]+1]=set(map(str_to_tuple,smart_split(row[1])))
  for i in range(max_deg):
    dels[i+1] |= dels[i]

  return dels

def get_Er_term(F,E,B, coe, r, B_gens, F_gens):
  conn = sqlite3.connect("cohomology.db", check_same_thread=False)
  conn.row_factory = sqlite3.Row
  cursor = conn.cursor()
  cursor.execute("SELECT r,kill FROM fibration WHERE F=? AND E=? AND B=? AND coe=?", (F,E,B, coe))
  result = cursor.fetchone()
  conn.commit()
  conn.close()

  max_p, max_q = 20, 20  # グリッドサイズ
  result_grid = [[[] for _ in range(max_q)] for _ in range(max_p)]

  # 1 ⊗ 1 = 1 を設定
  result_grid[0][0].append("1")

  dels = get_deleted(F,E,B,coe)

  for p in range(max_p):
    for q in range(max_q):
      if (p,q) in dels[r]: continue
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


def get_reference(F, E, B, coefficient):
  reference_row = df_reference[(df_reference["F"] == F) &
                                (df_reference["E"] == E) &
                                (df_reference["B"] == B) &
                                (df_reference["coe"] == int(coefficient))]
  return reference_row["reference"].values[0] if not reference_row.empty else None


# ホームページ (ファイブレーション選択)
@app.route("/", methods=["GET", "POST"])
def index():
  fibrations = get_fibrations()
  selected_fibration = None
  selected_coefficient = "1"
  cohomologies = None

  F,E,B="SU(3)","SU(4)","S^{7}" # 初期値
  cohomologies = get_fibration_cohomology((F, E, B), int(selected_coefficient))
  B_gens = get_cohomology_structure(B,selected_coefficient,"B")
  F_gens = get_cohomology_structure(F,selected_coefficient,"F")


  r=request.form.get("r", "2")
  tensor_product_grid = [[[] for _ in range(20)] for _ in range(20)]

  r=int(r)

  if request.method == "POST":
    selected_fibration = request.form.get("fibration")
    selected_coefficient = request.form.get("coefficient")
    
    if selected_fibration:
      F,E,B = smart_split(selected_fibration)
      cohomologies = get_fibration_cohomology((F, E, B), int(selected_coefficient))
      B_gens = get_cohomology_structure(B,selected_coefficient,"B")
      F_gens = get_cohomology_structure(F,selected_coefficient,"F")

      tensor_product_grid = get_Er_term(F,E,B, selected_coefficient, r, B_gens, F_gens)

  # print(get_ideal(F,selected_coefficient,"F"))
  # print(get_ideal(E,selected_coefficient,"E"))
  # print(get_ideal(B,selected_coefficient,"B"))

  E_gens = get_cohomology_structure(E, selected_coefficient, "E")
  reference = get_reference(F,E,B,selected_coefficient)

  return render_template("index.html", fibrations=fibrations \
    , selected_fibration=selected_fibration \
    , selected_coefficient=selected_coefficient \
    , cohomologies=cohomologies \
    , tensor_product_grid=tensor_product_grid \
    , E_gens=E_gens \
    , reference=reference \
    , r=r)

if __name__ == "__main__":
  app.run(debug=True)
