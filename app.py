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
df_differential = pd.read_csv("differential.csv").dropna()
df_differential2 = pd.read_csv("differential2.csv").dropna()

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

# differentialテーブル作成
cursor.execute("""
CREATE TABLE IF NOT EXISTS differential (
  F TEXT,
  E TEXT,
  B TEXT,
  coe INTEGER,
  r INTEGER,
  gen TEXT,
  dgen TEXT
)
""")
df_differential.to_sql("differential", conn, if_exists="replace", index=False)


# differential2 テーブル作成
cursor.execute("""
CREATE TABLE IF NOT EXISTS differential2 (
  F TEXT,
  E TEXT,
  B TEXT,
  coe INTEGER,
  r INTEGER,
  gen1 TEXT,
  gen2 TEXT,
  dgen1 TEXT,
  dgen2 TEXT
)
""")
df_differential2.to_sql("differential2", conn, if_exists="replace", index=False)

conn.commit()
conn.close()

print("CSVデータをSQLiteにインポートしました。")

max_deg = 20


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
        groups.append("\mathbb{Z}\{" + gen + "\}")
      else:
        groups.append(f"\mathbb{{Z}}_{{{gen_order[gen]}}}" + "\{" + gen + "\}")      

    groups.append("\cdots")
    cohomology_tex = r" \oplus ".join(groups) if groups else "0"
    return coe_tex, "\mathbb{Z}\{1\}\\oplus " + cohomology_tex

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

def get_differential(F,E,B,coe):
  conn = sqlite3.connect("cohomology.db", check_same_thread=False)
  conn.row_factory = sqlite3.Row
  cursor = conn.cursor()
  cursor.execute("SELECT r,gen,dgen FROM differential WHERE F=? AND E=? AND B=? AND coe=?", (F,E,B, coe))
  result = cursor.fetchall()
  conn.commit()
  conn.close()

  dif = [defaultdict(lambda: "0") for _ in range(max_deg+1)]

  for row in result:
    dif[row["r"]][row["gen"]] = row["dgen"]

  return dif

def get_differential2(F,E,B,coe):
  conn = sqlite3.connect("cohomology.db", check_same_thread=False)
  conn.row_factory = sqlite3.Row
  cursor = conn.cursor()
  cursor.execute("SELECT r,gen1,gen2,dgen1,dgen2 FROM differential2 WHERE F=? AND E=? AND B=? AND coe=?", (F,E,B, coe))
  result = cursor.fetchall()
  conn.commit()
  conn.close()

  dif = [defaultdict(lambda: "0") for _ in range(max_deg+1)]

  for row in result:
    dif[row["r"]][row["gen"]] = row["dgen"]

  return dif


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
    dels[row["r"]+1]=set(map(str_to_tuple,smart_split(row["kill"])))
  for i in range(max_deg):
    dels[i+1] |= dels[i]

  return dels

def get_target(B_gens, F_gens):
  target_grid = [[set() for _ in range(max_deg+1)] for _ in range(max_deg+1)]
  target_grid[0][0].add('1')

  for p in range(max_deg+1):
    for q in range(max_deg+1):
      B_list = B_gens[p] if p < len(B_gens) else []
      F_list = F_gens[q] if q < len(F_gens) else []
      for b in B_list:
        for f in F_list:
          if b == "1" and f == "1":
            continue  # すでに 1 を設定済み
          elif b == "1":
            target_grid[p][q].add(f"{f}")  # 1 ⊗ f = f
          elif f == "1":
            target_grid[p][q].add(f"{b}")  # b ⊗ 1 = b
          else:
            target_grid[p][q].add(f"{b} {f}")  # 一般形 b ⊗ f

  return target_grid

def extract_gen_and_exp(s):
  # 正規表現パターン：基本部分 (base) とオプションの上付き指数部分 (exponent)
  pattern = r"^(.*?)(?:\^\{(.*?)\})?$"
  match = re.match(pattern, s)
  if match:
    base = match.group(1)
    exp = match.group(2)
    if exp is None:
      return (base, 1)
    else:
      try:
        return (base, int(exp))
      except ValueError:
        return (base, exp)
  return (s, 1)

def str_to_exponent_map(x):
  res = defaultdict(int)
  split_x = x.split(' ')
  exps = list(map(extract_gen_and_exp,split_x))
  for s,e in exps:
    res[s] = e
  return res

def latex_product(x,y):
  if x == '0' or y == '0':
    return '0'
  if x == '1':
    if y == '1':
      return '1'
    else:
      return y
  else:
    if y == '1':
      return x
    else:
      x_exp = str_to_exponent_map(x)
      y_exp = str_to_exponent_map(y)
      xy_exp = defaultdict(int)
      for xi,exi in x_exp.items():
        xy_exp[xi] += exi
      for yi,eyi in y_exp.items():
        xy_exp[yi] += eyi
      xy_key = list(xy_exp.keys())
      xy_key.sort()
      res = []
      for xyi in xy_key:
        exyi = xy_exp[xyi]
        if exyi == 1:
          res.append(xyi)
        else:
          res.append(f"{xyi}^{{{exyi}}}")
      return ' '.join(res)

def get_tensor_tex(x):
  c,b,f=x
  if c==0 or b=='0' or f=='0':
    return '0'
  if c==1:
    return f"{b} \otimes {f}"
  return f"{c} ({b} \otimes {f})"

def latex_sum(x,y):
  c1,b1,f1=x
  c2,b2,f2=y
  if c1==0 or b1=='0' or f1=='0':
    if c2==0 or b2=='0' or f2=='0':
      return (0,'0','0')
    else:
      return y
  else:
    if c2==0 or b2=='0' or f2=='0':
      return x
    elif b1==b2 and f1==f2:
      return (c1+c2,b1,f1)
    else:
      return x


def get_Er_term(F,E,B, coe, r, B_gens, F_gens):
# SQLiteデータベースに接続し、微分情報を取得
  conn = sqlite3.connect("cohomology.db", check_same_thread=False)
  conn.row_factory = sqlite3.Row
  cursor = conn.cursor()
  cursor.execute("SELECT r,gen1,gen2,dgen1,dgen2 FROM differential2 WHERE F=? AND E=? AND B=? AND coe= ?", (F,E,B, coe))
  result = cursor.fetchall()
  conn.commit()
  conn.close()

# 各生成元の次数 (p, q) を格納する辞書
  degree = defaultdict(lambda: (0,0))
  
# 非ゼロ要素を格納するリストとセット
  non_zero_list = [[] for _ in range(r+1)] 
  non_zero_set = [set() for _ in range(r+1)]

# E_2-termの初期設定（B_gens と F_gens の直積を考える）
  for p in range(max_deg):
    for q in range(max_deg):
      for b in B_gens[p]:
        for f in F_gens[q]:
          degree[(b,f)] = (p,q)
          non_zero_list[2].append((b,f))

# スペクトル系列の結果を格納する3次元リスト
  result_grid = [[[[] for _ in range(max_deg+1)] \
                  for _ in range(max_deg+1)] for _ in range(max_deg)]

# E_2 のグリッドを設定
  for b,f in non_zero_list[2]:
    p,q = degree[(b,f)]
    result_grid[2][p][q].append(f"{b} \\otimes {f}")

# r = 2 の場合、E_2-term をそのまま返す
  if r == 2:
    return result_grid[2]

# 消滅する要素を管理するセットと微分辞書を初期化
  deleted_set = set()
  dif = [defaultdict(lambda: (0,'0','0')) for _ in range(r)]

# 取得したデータから微分を辞書に格納
  for row in result:
    if row["r"] < r:
      dif[row["r"]][(str(row["gen1"]),str(row["gen2"]))] \
        = (1,str(row["dgen1"]),str(row["dgen2"]))

# E_r-term を計算
  for i in range(2,r):
    non_zero_set[i] = set(non_zero_list[i])
# 微分による削除処理
    for (b1,f1) in non_zero_list[i]:
      p1,q1 = degree[(b1,f1)]
      for (b2,f2) in non_zero_list[i]:
        p2,q2 = degree[(b2,f2)]
        if p1+p2+i >= max_deg:
          continue  # 次数が範囲外の場合はスキップ

# d_r の積を計算
        c1,db1,df1 = dif[i][(b1,f1)]
        c2,db2,df2 = dif[i][(b2,f2)]
        b12 = latex_product(b1,b2)
        f12 = latex_product(f1,f2)
        if (b12,f12) not in non_zero_set[i]: continue
        (c12,bd12,fd12) = latex_sum((c1,latex_product(db1,b2),latex_product(df1,f2)), ((-1)**(p1+q1)*c2,latex_product(b1,db2),latex_product(f1,df2)))
        if coe!='0' and coe!='1':
          c12 %= int(coe)
        if coe=='1' and c12==-1:
          c12=1

# d_r 微分の結果が非ゼロなら削除対象に追加
        if c12!=0 and (bd12,fd12) in non_zero_set[i] and (bd12,fd12) != ('0','0'):
          dif[i][(b12,f12)] = (c12,bd12,fd12)
          deleted_set.add((b12,f12))
          if coe!='1' or (coe=='1' and c12==1):
            deleted_set.add((bd12,fd12))
          # if "^" in f12:
          print((b12,f12), c12,(bd12,fd12))
      # print(f"deleted_set = {deleted_set}")
# 削除されなかった要素を E_(r+1)-term に格納
    for b,f in non_zero_list[i]:
      if (b,f) in deleted_set:
        continue
      p,q = degree[(b,f)]
      result_grid[i+1][p][q].append(f"{b} \\otimes {f}")
      non_zero_list[i+1].append((b,f))

  return result_grid[r]


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
  tensor_product_grid = [[[] for _ in range(max_deg+1)] for _ in range(max_deg+1)]

  r=int(r)
  if r<2:
    r=2
  if r>max_deg-1:
    r=max_deg-1

  if request.method == "POST":
    selected_fibration = request.form.get("fibration")
    selected_coefficient = request.form.get("coefficient")
    
    if selected_fibration:
      F,E,B = smart_split(selected_fibration)
      cohomologies = get_fibration_cohomology((F, E, B), int(selected_coefficient))
      B_gens = get_cohomology_structure(B,selected_coefficient,"B")
      F_gens = get_cohomology_structure(F,selected_coefficient,"F")

      tensor_product_grid = get_Er_term(F,E,B,selected_coefficient,r,B_gens,F_gens)

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
  app.run(host="0.0.0.0", port=5000, debug=False)

