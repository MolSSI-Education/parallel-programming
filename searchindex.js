Search.setIndex({"docnames": ["00-distributed-memory", "00-shared-memory", "01-introduction", "02-distribution-intro", "03-distributed-examples-mpi4py", "04-distributed-examples", "05-shared-intro", "06-shared-examples", "index", "setup"], "filenames": ["00-distributed-memory.rst", "00-shared-memory.rst", "01-introduction.md", "02-distribution-intro.md", "03-distributed-examples-mpi4py.md", "04-distributed-examples.md", "05-shared-intro.md", "06-shared-examples.md", "index.rst", "setup.md"], "titles": ["Distributed Memory Parallelization", "Shared Memory Parallelization", "Introduction to Parallelization", "Introduction to Distributed-Memory Parallelization", "MPI Hands-On - mpi4py", "MPI Hands-On - C++", "Introduction to Shared-Memory Parallelization", "OpenMP Hands-On", "Parallel Programming", "Set Up"], "terms": {"introduct": [0, 1], "mpi": [0, 3, 8, 9], "hand": [0, 1, 2, 6, 8], "On": [0, 1, 8], "mpi4pi": [0, 8, 9], "c": [0, 2, 3, 4, 7, 8, 9], "openmp": [1, 6, 8, 9], "question": [2, 3, 4, 5, 6, 7, 8], "how": [2, 3, 4, 5, 7, 8], "doe": [2, 4, 5, 7, 8], "work": [2, 3, 4, 5, 7, 8, 9], "object": [2, 3, 4, 5, 6, 7, 8], "understand": [2, 3, 6, 7, 8], "motiv": [2, 8], "code": [2, 3, 4, 5, 7, 8], "machin": [2, 4, 5, 8, 9], "affect": [2, 8], "abil": [2, 8], "Be": [2, 8], "awar": [2, 8], "common": [2, 6, 8], "comput": [2, 3, 4, 5, 7, 8, 9], "chemistri": [2, 3, 8], "At": [2, 4, 5, 7], "some": [2, 4, 5, 6, 7, 9], "your": [2, 3, 4, 5, 7, 9], "career": 2, "you": [2, 3, 4, 5, 6, 7, 8, 9], "ve": [2, 7], "probabl": [2, 7], "ask": 2, "can": [2, 3, 4, 5, 6, 7, 8, 9], "make": [2, 4, 5, 7, 9], "my": [2, 8], "run": [2, 3, 4, 5, 6, 7, 8, 9], "faster": [2, 4, 5, 7], "Of": [2, 4, 5], "cours": [2, 4, 5, 8], "answer": 2, "thi": [2, 3, 4, 5, 7, 8, 9], "depend": [2, 6, 7, 9], "sensit": 2, "specif": [2, 3, 5, 6, 9], "situat": [2, 3], "here": [2, 4, 7, 8, 9], "ar": [2, 3, 4, 5, 6, 7, 8, 9], "few": [2, 4, 5, 7, 9], "thing": [2, 7], "might": [2, 3, 4, 5, 6, 7], "try": [2, 4, 5, 6, 7], "do": [2, 3, 4, 5, 6, 7, 8, 9], "optim": [2, 3, 7, 8], "move": [2, 7], "computation": 2, "demand": 2, "part": [2, 3, 4, 5, 7], "from": [2, 3, 4, 5, 6, 7, 9], "an": [2, 3, 4, 5, 6, 7, 8, 9], "interpret": 2, "languag": 2, "python": [2, 3, 4, 8, 9], "rubi": 2, "etc": [2, 3, 4, 5], "compil": [2, 5, 7, 8, 9], "julia": 2, "rust": 2, "us": [2, 3, 4, 5, 6, 7, 8, 9], "better": [2, 5, 7], "theoret": 2, "method": [2, 4, 5, 8], "requir": [2, 3, 4, 5, 6, 7, 8], "less": [2, 6], "same": [2, 3, 4, 5, 6, 7], "accuraci": 2, "each": [2, 3, 4, 5, 6, 7], "abov": [2, 3, 4, 5], "approach": [2, 3, 5, 6, 7], "intend": 2, "reduc": [2, 3, 7], "total": [2, 3, 4, 5, 7], "amount": [2, 3, 4, 5, 7], "A": [2, 3, 7], "differ": [2, 3, 4, 5, 6, 7], "strategi": [2, 3], "speed": [2, 4, 5, 6, 7], "up": [2, 4, 5, 7], "which": [2, 3, 4, 5, 6, 7, 8, 9], "split": 2, "among": [2, 3], "multipl": [2, 3, 4, 5, 6, 7], "process": [2, 3, 4, 5, 6], "unit": 2, "labor": 2, "simultan": [2, 4, 5, 6, 7], "The": [2, 3, 4, 5, 6, 7, 8, 9], "includ": [2, 3, 4, 5, 6, 7, 8], "central": 2, "cpu": 2, "graphic": 2, "gpu": 2, "vector": [2, 3, 4, 5, 7], "vpu": 2, "someth": [2, 4, 5, 7], "similar": [2, 7, 9], "just": [2, 4, 5, 7], "construct": 2, "worker": 2, "build": [2, 5, 7, 9], "hous": 2, "than": [2, 3, 4, 5, 6], "singl": [2, 4, 5, 6, 7], "complet": [2, 3, 4, 5, 6, 9], "calcul": [2, 3, 4, 5, 6, 7], "For": [2, 3, 4, 5, 6, 9], "exampl": [2, 3, 6], "take": [2, 4, 5, 7, 9], "1": [2, 3, 6, 9], "hour": 2, "one": [2, 4, 5, 6, 7], "possibl": [2, 4, 5], "across": [2, 4, 5, 7], "two": [2, 3, 4, 5, 6, 7], "onli": [2, 3, 4, 5, 6, 7], "30": [2, 7], "minut": [2, 9], "note": [2, 4, 5, 7], "cannot": 2, "fact": [2, 7], "gener": [2, 4, 5, 7], "introduc": [2, 3, 7], "addit": [2, 3, 4, 5], "associ": [2, 5], "commun": [2, 6, 8], "coordin": [2, 3, 4, 5, 7], "between": [2, 3, 4, 5], "In": [2, 3, 4, 5, 6, 7], "t": [2, 3, 4, 5, 7, 9], "serial": [2, 3, 5], "least": [2, 4, 5], "n": [2, 4, 5, 6, 7, 9], "principl": [2, 7], "more": [2, 3, 4, 5, 6, 7], "formal": 2, "express": 2, "through": [2, 4, 5, 7], "amdahl": 2, "": [2, 4, 5, 7, 8, 9], "law": 2, "primari": [2, 3, 4, 5], "goal": 2, "ensur": [2, 4, 5, 7], "actual": [2, 4, 5], "runtim": [2, 3], "close": [2, 7, 9], "ideal": [2, 7], "field": 2, "high": [2, 7], "perform": [2, 3, 4, 5, 7, 8], "concept": 2, "its": [2, 3, 4, 5, 6, 7], "logic": [2, 4, 5], "limit": [2, 3, 6, 7, 8], "rather": 2, "over": [2, 4, 5, 6, 7], "those": [2, 5], "find": [2, 7], "typic": [2, 4, 5], "desktop": 2, "laptop": [2, 9], "applic": [2, 5], "involv": [2, 3, 5, 6], "supercomput": 2, "consist": [2, 4, 5, 7], "mani": [2, 3, 4, 5, 7], "thousand": 2, "dizzi": 2, "scale": 2, "often": [2, 3, 4, 5, 7], "veri": [2, 6], "difficult": [2, 6], "although": [2, 3, 5, 6], "10": [2, 4, 7], "abl": [2, 6, 9], "time": [2, 4, 5, 6, 7, 9], "1000": [2, 4, 5], "most": [2, 3, 4, 5, 7], "them": [2, 4], "would": [2, 5, 7, 8, 9], "stand": [2, 4, 5], "around": [2, 7], "wait": 2, "unless": 2, "thei": [2, 4, 5, 6], "big": 2, "have": [2, 3, 4, 5, 6, 7, 9], "good": [2, 4, 5, 7, 9], "manag": 2, "whether": [2, 5, 6], "molecular": [2, 8], "properti": 2, "larg": [2, 3, 5, 7], "number": [2, 3, 4, 5, 6, 7], "quickli": [2, 7], "becom": 2, "overwhelm": 2, "purpos": 2, "lesson": [2, 3, 4, 5], "basic": [2, 7, 8], "behind": 2, "call": [2, 3, 4, 5, 7, 9], "attent": 2, "essenti": [2, 9], "help": [2, 4, 5, 7], "handl": [2, 3, 7, 8], "sometim": 2, "chaotic": 2, "confus": 2, "problem": [2, 3, 4, 5, 7], "softwar": [2, 8, 9], "first": [2, 3, 4, 5, 6, 7, 9], "level": [2, 7], "design": [2, 4, 5], "physic": [2, 3, 6], "imag": [2, 4], "below": [2, 3, 5, 9], "show": [2, 3], "rough": 2, "element": [2, 3, 4, 5, 7], "node": [2, 3, 6], "four": [2, 7], "shown": 2, "pictur": 2, "modern": 2, "interconnect": 2, "cluster": 2, "respect": [2, 7], "think": [2, 7], "independ": [2, 4, 5, 7], "connect": 2, "other": [2, 3, 4, 5, 6, 7], "within": [2, 3, 4, 5], "local": [2, 5, 6], "network": 2, "ha": [2, 3, 4, 5, 6, 7], "group": [2, 3], "core": [2, 3, 5, 6], "microprocessor": 2, "respons": [2, 7], "earlier": [2, 5], "dai": [2, 6], "all": [2, 3, 4, 5, 6, 7], "had": [2, 5], "todai": 2, "nearli": [2, 7], "multi": [2, 7], "share": [2, 7], "access": [2, 5, 6], "memori": [2, 7], "ram": [2, 6], "directli": [2, 6], "instead": [2, 3, 4, 5, 7], "three": [2, 7], "cach": [2, 6, 7], "e": [2, 4, 5, 6, 7], "g": 2, "l1": 2, "l2": 2, "exist": 2, "both": [2, 3, 6], "term": 2, "latenc": 2, "bandwidth": [2, 5], "smaller": [2, 3, 5], "capac": 2, "being": [2, 3, 4, 5, 7], "smallest": 2, "fastest": 2, "follow": [2, 3, 4, 5, 7, 9], "l3": 2, "when": [2, 4, 5, 6, 7, 9], "need": [2, 3, 4, 5, 6, 7, 9], "inform": [2, 3, 4, 5, 7], "search": 2, "If": [2, 3, 4, 5, 6, 7, 8, 9], "desir": 2, "found": [2, 9], "copi": [2, 3, 4, 5, 6, 7, 9], "ani": [2, 3, 4, 5, 6, 7, 9], "detail": [2, 3, 5, 6], "complex": [2, 3, 6], "won": [2, 3, 7], "expert": 2, "nonetheless": 2, "worthwhil": 2, "rate": 2, "import": [2, 3, 4, 5, 7, 9], "factor": [2, 6, 7], "effici": [2, 4, 5], "improv": [2, 3, 4, 5, 7], "easi": [2, 7], "hyper": 2, "focus": 2, "oper": [2, 3, 4, 5, 6, 7, 8, 9], "while": [2, 4, 5, 7], "ignor": [2, 5, 7], "distribut": [2, 6, 9], "case": [2, 3, 4, 5, 6, 7], "chang": [2, 4, 5, 7, 9], "alloc": [2, 3, 4, 5], "One": [2, 5, 6, 7, 9], "tool": 2, "assess": 2, "degre": 2, "bound": 2, "mathemat": [2, 4, 5], "versu": 2, "rooflin": 2, "model": 2, "also": [2, 3, 4, 5, 7, 9], "intra": 2, "fundament": [2, 8], "inter": [2, 6], "send": [2, 4, 5], "receiv": [2, 4, 5], "anoth": [2, 3, 4, 5, 7], "much": [2, 4, 5, 7], "It": [2, 3, 4, 5, 6, 7, 8], "increasingli": 2, "emploi": 2, "entir": [2, 4, 5, 6], "techniqu": [2, 3, 8], "instanc": [2, 3, 4, 5], "execut": [2, 3, 4, 5, 7], "own": [2, 3, 4, 5, 6, 7], "simul": [2, 4, 5, 7], "word": 2, "2": [2, 3, 6], "twice": 2, "mitig": 2, "thread": [2, 6, 7], "advantag": [2, 3, 6, 8], "becaus": [2, 3, 4, 5, 7], "must": [2, 3, 5, 6, 7], "support": 2, "instruct": [2, 5, 9], "data": [2, 4, 5, 6, 7], "simd": 2, "mean": [2, 4, 5], "operand": 2, "want": [2, 4, 5, 7], "multipli": [2, 3], "scalar": 2, "enabl": [2, 4, 5, 6], "individu": [2, 3, 4, 5, 6], "variou": 2, "wai": [2, 3, 4, 5, 7], "influenc": 2, "cover": 2, "heterogen": 2, "fpga": 2, "trend": 2, "recent": 2, "year": 2, "acceler": 2, "form": 2, "mutual": 2, "exclus": 2, "benefit": [2, 3], "notabl": 2, "popular": [2, 6], "These": [2, 5, 6, 7], "tutori": [2, 6], "focu": [2, 4, 5], "next": [2, 3, 4, 5, 6, 7, 9], "section": [2, 4, 5, 6, 7, 9], "we": [2, 3, 4, 5, 6, 7, 9], "subject": [2, 8], "activ": [2, 6, 9], "episod": [2, 3, 6, 8], "3": [2, 9], "4": [2, 4, 5, 8, 9], "5": [2, 4, 5, 7], "6": [2, 5, 7, 8], "link": [2, 3, 5, 6, 7, 9], "md": [2, 3, 5, 6, 7], "achiev": [2, 4, 5, 7], "There": [2, 4, 5], "hardwar": 2, "avail": 2, "what": [3, 4, 5, 6, 7, 8], "i": [3, 4, 5, 6, 7, 8, 9], "illustr": 3, "figur": 3, "identif": 3, "known": [3, 4, 5, 7], "rank": [3, 4, 5], "assign": [3, 4, 5], "sequenti": [3, 7], "increas": [3, 4, 5], "order": [3, 4, 5, 7], "start": [3, 6, 7, 9], "zero": [3, 4, 7], "exact": [3, 6, 7, 9], "duplic": [3, 4, 5], "write": [3, 6, 7, 8], "without": [3, 4, 5], "thought": [3, 5], "simpli": 3, "base": [3, 9], "suppos": [3, 6], "dot": 3, "product": [3, 4, 5], "b": [3, 4, 5, 7], "look": [3, 4, 5, 7, 9], "like": [3, 4, 5, 7, 8], "togeth": 3, "second": [3, 5, 6, 7], "add": [3, 4, 5, 6, 7], "previou": [3, 4], "result": [3, 4, 5, 6, 7], "so": [3, 4, 5, 7], "nproc": [3, 4, 5], "my_rank": [3, 4, 5], "could": [3, 4, 5, 6, 7], "th": 3, "issu": [3, 7, 8], "now": [3, 4, 5, 7, 9], "somehow": 3, "piec": 3, "get": [3, 6, 7], "accomplish": 3, "until": 3, "alreadi": [3, 4, 5, 7], "begin": [3, 4, 5, 7], "see": [3, 4, 5, 7, 8, 9], "challeng": [3, 6], "implement": [3, 4, 5, 6, 7, 8], "divid": 3, "stitch": 3, "back": [3, 6, 7], "third": [3, 7], "major": [3, 4, 5], "keep": [3, 4, 5], "usag": 3, "accept": [3, 4, 5], "tend": [3, 4, 5, 6], "store": [3, 5, 6, 7], "redund": 3, "nuclear": [3, 4, 5], "atom": [3, 4, 5, 7], "compris": 3, "system": [3, 5, 7, 9], "quantiti": 3, "By": [3, 5], "default": [3, 5, 6, 7], "uniqu": [3, 7], "modifi": [3, 4, 5, 7], "behavior": [3, 4, 5], "intellig": [3, 4, 5], "written": 3, "isn": [3, 4, 5, 7, 9], "necessari": [3, 4, 5], "entireti": 3, "arrai": [3, 4, 5, 7], "real": [3, 7], "world": [3, 7], "far": [3, 6], "easili": 3, "effort": [3, 4, 5], "commonli": 3, "messag": [3, 4, 5, 7, 9], "pass": [3, 4, 5], "interfac": [3, 4, 5], "demonstr": 3, "mai": [3, 7, 9], "whichev": [3, 4, 5], "relev": 3, "mechan": 3, "acheiv": 3, "larger": [3, 5, 7], "parallel": [4, 5, 7], "learn": [4, 5, 7, 8], "prepar": [4, 8], "environ": [4, 5, 6, 8], "explor": [4, 6, 8], "ll": [4, 5], "example1": [4, 7], "simpl": [4, 5, 7], "__name__": 4, "__main__": 4, "print": [4, 5, 7, 9], "acquir": [4, 5], "file": [4, 5, 7, 8, 9], "shell": [4, 5, 7, 9], "git": [4, 5], "clone": [4, 5], "github": [4, 5, 7], "com": [4, 5, 9], "molssi": [4, 5, 8], "educ": [4, 5, 8], "program": [4, 5], "cd": [4, 5, 7], "py": 4, "output": [4, 5, 7, 9], "let": [4, 5, 7, 9], "done": [4, 5, 9], "mpiexec": [4, 5, 9], "command": [4, 5, 7, 9], "provid": [4, 5], "mpirun": [4, 5], "usual": [4, 5, 7], "alwai": [4, 5, 7], "whenev": [4, 5], "should": [4, 5, 7, 8, 9], "guarante": [4, 5], "v": [4, 5, 7], "standard": [4, 5], "varieti": [4, 5], "architectur": [4, 5, 6, 8], "defin": [4, 5, 7], "syntax": [4, 5], "semant": [4, 5], "librari": [4, 5], "routin": [4, 5], "openmpi": [4, 5], "mpich": [4, 5, 9], "m": [4, 5, 9], "option": [4, 5, 9], "technic": [4, 5], "doesn": [4, 5, 7], "either": [4, 5], "describ": [4, 5], "guidelin": [4, 5], "prefer": [4, 5, 9], "format": [4, 5], "lanch": [4, 5], "number_of_process": [4, 5], "command_to_launch_cod": [4, 5], "launch": [4, 5], "long": [4, 5, 7, 9], "processor": [4, 5], "howev": [4, 5], "certain": [4, 5], "variabl": [4, 5, 6, 7], "argument": [4, 5], "obviou": [4, 5], "yet": [4, 5, 7], "aren": [4, 5], "unawar": [4, 5], "via": [4, 5, 7], "span": [4, 5], "comm_world": [4, 9], "class": 4, "queri": [4, 5], "about": [4, 5, 6, 7], "function": [4, 5, 7], "edit": [4, 5, 7, 9], "read": [4, 5, 7, 9], "world_comm": [4, 5], "world_siz": [4, 5], "get_siz": [4, 9], "get_rank": 4, "size": [4, 5, 7], "str": 4, "Then": [4, 5, 7], "tell": [4, 5, 7], "u": [4, 5, 7], "uniq": [4, 5], "integ": [4, 5, 6], "rang": [4, 5], "0": [4, 5, 6, 7], "allow": [4, 5, 6, 7], "identifi": [4, 5, 6, 7, 8], "return": [4, 5, 7], "go": [4, 5, 7], "ahead": [4, 5, 7], "As": [4, 5, 6, 7], "told": [4, 5, 7], "don": [4, 5, 7, 9], "necessarili": [4, 5], "out": [4, 5, 7], "reach": 4, "again": [4, 5, 7], "thier": [4, 5], "rerun": [4, 5], "valu": [4, 5, 6, 7], "script": [4, 7, 8], "example2": [4, 7], "math": [4, 5, 7], "numpi": [4, 9], "averag": [4, 5, 7], "5000001": 4, "account": [4, 5], "timer": [4, 5], "wtime": 4, "current": [4, 5], "walltim": [4, 5], "determin": [4, 5, 7], "spent": [4, 5, 7], "initi": [4, 5, 7], "start_tim": [4, 5, 7], "np": 4, "ones": 4, "end_tim": [4, 5], "indic": [4, 5], "realli": [4, 5, 7], "everi": [4, 5, 7], "sinc": [4, 5], "messi": [4, 5], "line": [4, 5, 7, 9], "top": [4, 5, 7], "intial": [4, 5], "final": [4, 5, 6, 7], "10000000": 4, "sum": [4, 5, 7], "03975701332092285": 4, "569957971572876": 4, "173098087310791": 4, "609341859817505": 4, "042365074157714844": 4, "9863519668579102": 4, "9583611488342285": 4, "9468209743499756": 4, "cooper": [4, 5], "decid": [4, 5, 7], "workload": [4, 5], "befor": [4, 5, 6, 7, 9], "my_start": [4, 5], "my_end": [4, 5], "repres": [4, 5], "updat": [4, 5, 6, 7], "loop": [4, 5, 7], "global": [4, 5], "To": [4, 5, 7, 8, 9], "replac": [4, 5, 6, 7], "world_sum": 4, "sum_np": 4, "empti": 4, "recv": [4, 5], "doubl": [4, 5, 7], "sourc": [4, 5], "tag": [4, 5], "77": [4, 5], "els": [4, 5], "dest": [4, 5], "paramet": [4, 5], "type": [4, 5, 7, 8, 9], "precis": [4, 5], "datatyp": [4, 5], "consult": [4, 5], "byte": 4, "8": [4, 5, 7], "binari": [4, 5, 9], "digit": [4, 5], "char": [4, 5], "unsigned_char": 4, "unsign": [4, 5], "short": [4, 5, 7], "sign": [4, 5], "int": [4, 5, 7], "unsigned_short": 4, "unsigned_long": 4, "float": [4, 5], "04637002944946289": 4, "9484930038452148": 4, "914314031600952": 4, "6889588832855225": 4, "inde": [4, 5], "gone": [4, 5], "down": [4, 5], "easier": [4, 5, 7], "iter": [4, 5, 7], "04810309410095215": 4, "0196259021759033": 4, "2053139209747314": 4, "721329927444458": 4, "nice": [4, 5, 7], "surprisingli": [4, 5], "enough": [4, 5, 7], "expens": [4, 5, 7], "calat": [4, 5], "04351997375488281": 4, "503791093826294": 4, "2048840522766113": 4, "7626049518585205": 4, "thank": [4, 5], "ad": [4, 5, 7], "care": [4, 5], "stop": [4, 5], "realiti": [4, 5], "though": [4, 5], "resourc": [4, 5], "concern": [4, 5], "even": [4, 5, 6, 7], "made": [4, 5, 6, 7], "decreas": [4, 5], "That": [4, 5, 7], "our": [4, 5, 7], "correspond": [4, 5], "reason": [4, 5, 6], "ever": [4, 5], "009948015213012695": 4, "5988950729370117": 4, "2081310749053955": 4, "7307591438293457": 4, "previous": [4, 5], "autom": [4, 5], "complic": [4, 5], "particular": [4, 5], "op": [4, 5], "root": [4, 5], "specifi": [4, 7], "set": [4, 5, 7], "caus": [4, 5, 7], "onto": [4, 5, 9], "descript": [4, 5], "max": [4, 5], "maximum": [4, 5], "min": [4, 5], "minimum": [4, 5], "prod": 4, "land": 4, "AND": [4, 5], "band": 4, "bit": [4, 5], "wise": [4, 5], "mpi_byt": [4, 5], "lor": 4, "OR": [4, 5], "bor": 4, "lxor": 4, "xor": [4, 5], "bxor": 4, "maxloc": 4, "locat": [4, 5, 6, 7], "minloc": 4, "simpler": [4, 5, 9], "view": [4, 5], "example3": [4, 7], "mont": [4, 5], "carlo": [4, 5], "248": 4, "52688099543923": 4, "2000": 4, "588491394826892": 4, "3000": 4, "9309007491547571": 4, "4000": 4, "8247648102916196": 4, "5000": 4, "715929587912762": 4, "6000": 4, "362217832200815": 4, "7000": 4, "570585267104749": 4, "8000": 4, "649439720181915": 4, "9000": 4, "65428738463388": 4, "10000": [4, 5], "73417919011543": 4, "21": [4, 7], "389078855514526": 4, "energi": [4, 5, 7], "013432502746582": 4, "decis": [4, 5], "09333038330078125": 4, "vast": [4, 5], "get_particle_energi": [4, 5], "where": [4, 5, 7], "def": 4, "box_length": [4, 5], "i_particl": [4, 5], "cutoff2": [4, 5], "pairwis": 4, "lennard": [4, 5], "jone": [4, 5], "particl": [4, 5, 7], "period": 4, "box": 4, "r_i": 4, "list": [4, 7], "potit": 4, "vection": 4, "r_j": 4, "j": [4, 7], "length": [4, 7], "rij2": [4, 5], "squar": 4, "shortest": 4, "distanc": 4, "e_tot": [4, 5], "i_posit": [4, 5], "particle_count": [4, 5], "len": 4, "j_particl": [4, 5], "j_posit": [4, 5], "minimum_image_dist": [4, 5], "e_pair": 4, "lennard_jones_potenti": [4, 5], "fairli": [4, 5], "straightforward": [4, 5, 7], "interact": [4, 5, 7], "pair": [4, 5, 7], "subset": [4, 5, 6], "know": [4, 5], "comm": [4, 5], "main": [4, 5, 7], "current_energi": [4, 5], "simulation_cutoff2": [4, 5], "place": [4, 5], "stride": [4, 5], "offset": [4, 5], "12": [4, 5, 7], "9": [4, 5, 7], "13": [4, 5], "e_singl": 4, "e_sum": [4, 5], "35480909996": 4, "566864": 4, "66252436255523": 4, "72": 4, "86936127660856": 4, "08": 4, "93141042416342": 4, "66": 4, "256171999678073": 4, "88": 4, "162015453630529e": 4, "1620181302289283e": 4, "162018130377518e": 4, "1620181324457333e": 4, "1620182854716e": 4, "31": 4, "748733043670654": 4, "112581253051758": 4, "21792912483215332": 4, "seem": [4, 5], "right": [4, 5], "went": [4, 5], "wrong": [4, 5], "none": [4, 5, 7], "allreduc": 4, "5402881": 4, "246788438": 4, "5403807": 4, "559181325": 4, "5403898": 4, "801044374": 4, "5403916": 4, "261693102": 4, "5403921": 4, "433225453": 4, "5403923": 4, "534017933": 4, "5403924": 4, "646963553": 4, "5403925": 4, "292483066": 4, "63053995": 4, "5403926": 4, "272461226": 4, "43": 4, "26621890068054": 4, "42": 4, "664116621017456": 4, "16298675537109375": 4, "still": [4, 5], "randomli": [4, 5], "displac": [4, 5], "select": [4, 5], "contribut": [4, 5, 7, 8], "configur": [4, 5], "end": [4, 5, 7], "fix": [4, 5, 6, 7, 8], "broadcast": [4, 5], "sync": [4, 5], "generate_initial_st": [4, 5], "build_method": 4, "num_particl": [4, 5], "bcast": 4, "i_step": [4, 5], "n_step": [4, 5], "n_trial": [4, 5], "random": [4, 5, 7], "randint": 4, "random_displac": [4, 5], "rand": 4, "max_displac": [4, 5], "i_particle_buf": 4, "start_decision_tim": [4, 5], "delta_": [4, 5], "proposed_energi": [4, 5], "accept_or_reject": [4, 5], "beta": [4, 5], "total_pair_energi": [4, 5], "n_accept": [4, 5], "total_energi": [4, 5], "tail_correct": [4, 5], "energy_arrai": [4, 5], "mod": 4, "freq": [4, 5], "tune_displac": [4, 5], "adjust_displac": [4, 5], "total_decision_tim": [4, 5], "52688099525105": 4, "588491394638726": 4, "9309007493429244": 4, "824764810479789": 4, "715929588100931": 4, "3622178323889855": 4, "570585267292914": 4, "649439720370088": 4, "65428738482205": 4, "734179190303595": 4, "671964883804321": 4, "892877340316772": 4, "15127253532409668": 4, "expect": [4, 5], "simplic": [4, 5], "choic": [4, 5], "proper": [5, 8], "non": [5, 8], "block": [5, 7, 8], "debugg": [5, 8], "cpp": [5, 7], "iostream": 5, "argc": 5, "argv": 5, "std": [5, 7], "cout": 5, "endl": 5, "mkdir": 5, "cmake": [5, 7, 9], "dcmake_c_compil": 5, "mpicc": 5, "dcmake_cxx_compil": 5, "mpicxx": 5, "dcmake_fortran_compil": 5, "mpifort": 5, "callout": 5, "mpi_comm_world": 5, "h": [5, 7], "mpi_init": 5, "mpi_comm_s": 5, "mpi_comm_rank": 5, "after": [5, 7, 9], "mpi_fin": 5, "recompil": 5, "header": 5, "forc": [5, 7], "exit": 5, "never": [5, 6], "termin": [5, 8, 9], "leav": 5, "anyth": [5, 7], "indefint": 5, "re": 5, "wast": [5, 7], "massiv": 5, "fail": 5, "mpi_abort": 5, "goe": 5, "mind": 5, "succesfulli": 5, "equal": [5, 7], "mpi_success": 5, "wa": [5, 7], "successfulli": 5, "otherwis": [5, 7], "automat": 5, "abort": 5, "encount": [5, 7], "yourself": 5, "mpi_errhandler_set": 5, "mpi_errors_return": 5, "check": [5, 7, 9], "100000001": 5, "mpi_wtim": 5, "new": [5, 7], "200000000": 5, "delet": 5, "544075": 5, "624939": 5, "258915": 5, "266418": 5, "640894": 5, "893775": 5, "38309": 5, "330192": 5, "wors": 5, "compet": 5, "seriou": 5, "due": 5, "extrem": [5, 7], "manipul": 5, "explicitli": 5, "reciev": 5, "mpi_send": 5, "mpi_recv": 5, "const": [5, 7], "void": [5, 7], "buf": 5, "count": [5, 7], "mpi_datatyp": 5, "mpi_comm": 5, "pointer": 5, "buffer": 5, "sent": 5, "destin": 5, "mpi_statu": 5, "statu": 5, "hold": 5, "mpi_any_sourc": 5, "match": 5, "mpi_any_tag": 5, "structur": 5, "ran": 5, "immedi": 5, "partial_averag": 5, "mpi_doubl": 5, "mpi_char": 5, "mpi_unsigned_char": 5, "mpi_short": 5, "mpi_unsigned_short": 5, "mpi_int": 5, "mpi_unsign": 5, "mpi_long": 5, "mpi_unsigned_long": 5, "mpi_float": 5, "mpi_pack": 5, "mpi_unpack": 5, "63251": 5, "31379": 5, "89099": 5, "100575": 5, "636685": 5, "66542": 5, "466888": 5, "0871116": 5, "initialz": 5, "159471": 5, "183946": 5, "193497": 5, "0847806": 5, "16013": 5, "176896": 5, "190774": 5, "0871552": 5, "particularli": [5, 6], "mpi_reduc": 5, "sendbuf": 5, "recvbuf": 5, "mpi_op": 5, "address": 5, "mpi_max": 5, "mpi_min": 5, "mpi_sum": 5, "mpi_prod": 5, "mpi_land": 5, "mpi_band": 5, "mpi_lor": 5, "mpi_bor": 5, "mpi_lxor": 5, "mpi_bxor": 5, "mpi_maxloc": 5, "mpi_minloc": 5, "produc": [5, 6, 7], "mc": 5, "497000": 5, "28643": 5, "498000": 5, "28989": 5, "499000": 5, "96743": 5, "500000": [5, 7], "06861": 5, "59121": 5, "47059": 5, "0425119": 5, "16000": 5, "19516e": 5, "19": 5, "17000": 5, "19517e": 5, "bad": 5, "OF": 5, "ONE": 5, "pid": 5, "81043": 5, "AT": 5, "taylor": 5, "macbook": 5, "pro": 5, "clean": 5, "remain": 5, "THE": 5, "cleanup": 5, "WITH": 5, "string": 5, "segment": 5, "fault": 5, "11": 5, "signal": 5, "refer": 5, "pleas": [5, 8, 9], "faq": 5, "page": 5, "debug": [5, 6], "suggest": 5, "mpi_allreduc": 5, "28644": 5, "2899": 5, "96744": 5, "06862": 5, "38658": 5, "24201": 5, "0563176": 5, "certainli": 5, "experi": 5, "Near": 5, "seed": 5, "pre": 5, "mt19937": 5, "mt": 5, "random_devic": 5, "rd": 5, "switch": [5, 7], "7": [5, 7, 9], "9279e": 5, "73895": 5, "59179": 5, "0549728": 5, "crazi": 5, "unphys": 5, "confirm": 5, "give": [5, 7], "lead": [5, 6], "utter": 5, "chao": 5, "throughout": 5, "tempt": 5, "sound": [5, 6], "open": [5, 7, 8, 9], "diverg": 5, "rememb": 5, "infinit": 5, "accur": 5, "slight": 5, "descrep": 5, "happen": [5, 7], "given": [5, 9], "mpi_bcast": 5, "floor": 5, "dist": 5, "test": 5, "reject": 5, "step": [5, 7], "bool": 5, "revert": 5, "posit": 5, "round": 5, "06666": 5, "10058": 5, "98052": 5, "95301": 5, "16": [5, 7], "7881": 5, "15": [5, 7], "9948": 5, "0690625": 5, "consider": [5, 6], "were": 5, "small": 5, "under": [5, 8], "setup": 5, "comment": [5, 7], "100": [5, 7], "overhead": [5, 7], "extra": [5, 7], "somewhat": [5, 7], "100000": 5, "97000": 5, "612": 5, "067": 5, "98000": 5, "609": 5, "113": 5, "99000": 5, "603": 5, "538": 5, "599": 5, "461": 5, "41": 5, "0191": 5, "39": 5, "9748": 5, "011933": 5, "99": [5, 7], "1126": 5, "93": 5, "454": 5, "91": 5, "1246": 5, "87": 5, "397": 5, "22": 5, "2873": 5, "4661": 5, "0175401": 5, "clearli": [5, 7], "speedup": 5, "spawn": 6, "privat": 6, "worri": 6, "conveni": 6, "bizarr": 6, "bug": 6, "clear": [6, 7], "mention": 6, "fetch": 6, "moment": [6, 7], "subsequ": 6, "But": 6, "numer": 6, "uncontrol": 6, "programm": 6, "algorithm": 6, "microscop": 6, "defect": 6, "minor": 6, "regard": 6, "liter": 6, "hot": 6, "cold": 6, "maintain": 6, "assum": [6, 8], "despit": 6, "subtl": 6, "signific": [6, 7], "wherea": 6, "regardless": 6, "dure": 6, "lower": 6, "race": [7, 8], "condit": [7, 8], "directori": 7, "repositori": [7, 8], "omp": 7, "text": [7, 9], "editor": [7, 9], "hello": 7, "stdio": 7, "printf": 7, "sh": [7, 9], "turn": 7, "direct": 7, "pragma": 7, "pragmat": 7, "put": 7, "curli": 7, "bracket": 7, "distinct": 7, "rebuild": 7, "special": 7, "omp_num_thread": 7, "export": [7, 9], "obvious": 7, "thread_id": 7, "omp_get_thread_num": 7, "finish": 7, "slightli": 7, "plai": 7, "insid": 7, "num_thread": 7, "goodby": 7, "master": 7, "fork": 7, "join": 7, "resum": 7, "creat": [7, 8], "http": [7, 9], "en": 7, "wikipedia": 7, "org": 7, "wiki": 7, "e2": 7, "80": 7, "93join_model": 7, "media": 7, "fork_join": 7, "svg": 7, "500000001": 7, "sort": 7, "costli": 7, "omp_get_wtim": 7, "f": 7, "816225": 7, "statement": [7, 8], "402934": 7, "077723": 7, "900899": 7, "381398": 7, "781963": 7, "id": 7, "nthread": 7, "istart": 7, "iend": 7, "nthr": 7, "omp_get_num_thread": 7, "348737": 7, "sai": [7, 9], "328928": 7, "397873": 7, "076062": 7, "342933": 7, "123558": 7, "156250000": 7, "375000": 7, "960761": 7, "strong": 7, "hint": 7, "come": 7, "reduct": 7, "claus": 7, "383079": 7, "074402": 7, "346857": 7, "123527": 7, "948352": 7, "coupl": 7, "bottlenect": 7, "391194": 7, "439494": 7, "296057": 7, "105363": 7, "252643": 7, "344843": 7, "350311": 7, "188084": 7, "085440": 7, "987741": 7, "notic": 7, "last": 7, "touch": 7, "999": 7, "92079": 7, "129718": 7, "pe": 7, "16253": 7, "127101": 7, "000707": 7, "calc": 7, "992926": 7, "veloc": 7, "014128": 7, "coord": 7, "001717": 7, "014781": 7, "region": 7, "start_loop": 7, "natom": 7, "dx": 7, "dy": 7, "dr2": 7, "dr": 7, "sqrt": 7, "fx": 7, "fy": 7, "potenti": 7, "inner": 7, "286004": 7, "735108": 7, "18173": 7, "082058": 7, "001293": 7, "14": 7, "219656": 7, "014499": 7, "005649": 7, "246469": 7, "unfortun": 7, "longer": 7, "why": 7, "innermost": 7, "subtract": 7, "attempt": 7, "subtrat": 7, "outcom": 7, "overwrit": 7, "scenario": 7, "sever": 7, "elimin": 7, "avoid": 7, "error": [7, 8], "correct": [7, 9], "substanti": 7, "000792": 7, "954821": 7, "003412": 7, "005082": 7, "969381": 7, "been": 7, "said": 7, "enter": 7, "later": 7, "cost": 7, "vari": 7, "microsecond": 7, "spend": 7, "approxim": 7, "decent": 7, "outer": 7, "000737": 7, "770037": 7, "002868": 7, "005819": 7, "784812": 7, "definit": 7, "appar": 7, "half": 7, "realist": 7, "example4": 7, "03": 7, "84389e": 7, "17": 7, "11479": 7, "10384": 7, "14405": 7, "16274": 7, "17473": 7, "18386": 7, "19046": 7, "19516": 7, "19881": 7, "20184": 7, "20445": 7, "ke": 7, "p": 7, "00": 7, "2937": 7, "72894": 7, "18197": 7, "02586": 7, "15259": 7, "29692": 7, "371": 7, "04086690": 7, "3426": 7, "85950": 7, "18452": 7, "47277": 7, "15025": 7, "61327": 7, "408": 7, "04039200": 7, "3603": 7, "35315": 7, "18693": 7, "26130": 7, "15089": 7, "90815": 7, "441": 7, "03943282": 7, "3615": 7, "37003": 7, "19040": 7, "76407": 7, "15425": 7, "39404": 7, "436": 7, "03848593": 7, "3686": 7, "47202": 7, "19111": 7, "90579": 7, "43377": 7, "453": 7, "03809054": 7, "18": [7, 9], "3702": 7, "65533": 7, "19128": 7, "06154": 7, "40621": 7, "465": 7, "03778621": 7, "3812": 7, "43640": 7, "19237": 7, "88713": 7, "45072": 7, "472": 7, "03758076": 7, "24": 7, "3849": 7, "80094": 7, "19275": 7, "24791": 7, "44697": 7, "480": 7, "03741216": 7, "27": 7, "3962": 7, "65623": 7, "19388": 7, "13553": 7, "47930": 7, "492": 7, "03713681": 7, "3973": 7, "40391": 7, "19398": 7, "90927": 7, "50536": 7, "495": 7, "03703538": 7, "85": 7, "neigh": 7, "29": 7, "32": 7, "primarili": 7, "neighbor_list": 7, "consid": 7, "neighbor": 7, "timestep": 7, "evalu": 7, "neight": 7, "const_iter": 7, "ij": 7, "littl": 7, "awkward": 7, "separ": 7, "neighborlist": 7, "thrneight": 7, "cc": 7, "typdef": 7, "typedef": 7, "contain": 7, "appropri": 7, "thrneigh": 7, "coordt": 7, "ithr": 7, "virial": 7, "dt": 7, "prev": 7, "1e99": 7, "600": 7, "nneigh": 7, "relax": 7, "guess": 7, "potential_energi": 7, "virial_step": 7, "40320": 7, "90856": 7, "50": 7, "06": 7, "67": 7, "readi": 7, "subroutin": 7, "increment": 7, "dfx": 7, "dfy": 7, "unlik": 7, "restructur": 7, "portion": 7, "manual": 7, "declar": 7, "f_thread": 7, "xyt": 7, "virial_thread": 7, "pe_thread": 7, "vij": 7, "critic": [7, 9], "40314": 7, "90851": 7, "26": 7, "23": 7, "onc": [7, 9], "convert": 7, "creation": 7, "65": 7, "49": 7, "58": 7, "load": 7, "balanc": 7, "barrier": 7, "target": 7, "per": 7, "npair": 7, "steal": 7, "break": 7, "push_back": 7, "pop_back": 7, "donat": 7, "59": 7, "37": 7, "lot": 7, "resiz": 7, "slow": 7, "pairt": 7, "reserv": 7, "best": 7, "55": 7, "54": 7, "carefulli": 7, "refactor": 7, "scienc": 8, "institut": 8, "teach": 8, "emphasi": 8, "familiar": 8, "full": 8, "mission": 8, "continu": 8, "develop": 8, "report": 8, "submit": 8, "pull": 8, "request": 8, "student": 8, "window": 8, "navig": 8, "bash": [8, 9], "titl": 8, "download": [8, 9], "memoeri": 8, "point": 8, "collect": 8, "strongli": 9, "recommend": 9, "subsystem": 9, "ubuntu": 9, "04": 9, "menu": 9, "choos": 9, "usernam": 9, "password": 9, "wget": 9, "repo": 9, "archiv": 9, "anaconda3": 9, "2020": 9, "02": 9, "x86_64": 9, "reopen": 9, "init": 9, "path": 9, "echo": 9, "home": 9, "your_usernam": 9, "bin": 9, "bashrc": 9, "verifi": 9, "jupyt": 9, "notebook": 9, "browser": 9, "url": 9, "outlin": 9, "red": 9, "past": 9, "sudo": 9, "apt": 9, "vscode": 9, "extens": 9, "click": 9, "websit": 9, "sure": 9, "version": 9, "user": 9, "maco": 9, "xcode": 9, "everyon": 9, "great": 9, "deal": 9, "isol": 9, "stack": 9, "rest": 9, "molssi_pp": 9, "forg": 9, "everyth": 9, "openrt": 9, "packag": 9, "correctli": 9}, "objects": {}, "objtypes": {}, "objnames": {}, "titleterms": {"distribut": [0, 3, 8], "memori": [0, 1, 3, 4, 5, 6, 8], "parallel": [0, 1, 2, 3, 6, 8, 9], "share": [1, 6, 8], "introduct": [2, 3, 6, 8], "overview": [2, 3, 4, 5, 6, 7], "what": 2, "i": 2, "hpc": 2, "architectur": 2, "type": 2, "kei": [2, 3, 4, 5, 6, 7], "point": [2, 3, 4, 5, 6, 7], "mpi": [4, 5], "hand": [4, 5, 7], "On": [4, 5, 7], "mpi4pi": 4, "1": [4, 5, 7], "exampl": [4, 5, 7], "write": [4, 5], "hello": [4, 5], "world": [4, 5], "get": [4, 5], "start": [4, 5], "2": [4, 5, 7], "basic": [4, 5], "infrastructur": [4, 5], "commun": [4, 5], "reduc": [4, 5], "footprint": [4, 5], "collect": [4, 5], "3": [4, 5, 7], "c": 5, "error": 5, "handl": 5, "openmp": 7, "4": 7, "program": [8, 9], "prerequisit": 8, "workshop": 8, "lesson": 8, "set": [8, 9], "up": [8, 9], "anaconda": 9, "instal": 9, "window": 9, "wsl": 9, "without": 9, "mac": 9, "o": 9, "linux": 9, "creat": 9, "conda": 9, "environ": 9, "confirm": 9}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 8, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.viewcode": 1, "sphinx.ext.intersphinx": 1, "sphinx": 57}, "alltitles": {"Distributed Memory Parallelization": [[0, "distributed-memory-parallelization"]], "Shared Memory Parallelization": [[1, "shared-memory-parallelization"]], "Introduction to Parallelization": [[2, "introduction-to-parallelization"]], "Overview": [[2, null], [3, null], [4, null], [5, null], [6, null], [7, null]], "What is Parallelization?": [[2, "what-is-parallelization"]], "HPC Architecture": [[2, "hpc-architecture"]], "Types of Parallelization": [[2, "types-of-parallelization"]], "Key Points": [[2, null], [3, null], [4, null], [5, null], [6, null], [7, null]], "Introduction to Distributed-Memory Parallelization": [[3, "introduction-to-distributed-memory-parallelization"]], "MPI Hands-On - mpi4py": [[4, "mpi-hands-on-mpi4py"]], "1. Example 1": [[4, "example-1"], [5, "example-1"]], "Writing Hello World": [[4, "writing-hello-world"], [5, "writing-hello-world"]], "Getting Started with MPI": [[4, "getting-started-with-mpi"], [5, "getting-started-with-mpi"]], "Example 2": [[4, "example-2"], [5, "example-2"], [7, "example-2"]], "Basic Infrastructure": [[4, "basic-infrastructure"], [5, "basic-infrastructure"]], "Point-to-Point Communication": [[4, "point-to-point-communication"], [5, "point-to-point-communication"]], "Reducing the Memory Footprint": [[4, "reducing-the-memory-footprint"], [5, "reducing-the-memory-footprint"]], "Collective Communication": [[4, "collective-communication"], [5, "collective-communication"]], "Example 3": [[4, "example-3"], [5, "example-3"], [7, "example-3"]], "MPI Hands-On - C++": [[5, "mpi-hands-on-c"]], "Error Handling with MPI": [[5, "error-handling-with-mpi"]], "Introduction to Shared-Memory Parallelization": [[6, "introduction-to-shared-memory-parallelization"]], "OpenMP Hands-On": [[7, "openmp-hands-on"]], "Example 1": [[7, "example-1"]], "Example 4": [[7, "example-4"]], "Parallel Programming": [[8, "parallel-programming"]], "Prerequisites": [[8, null]], "Workshop Lessons": [[8, "workshop-lessons"]], "Set-Up": [[8, "set-up"]], "Introduction": [[8, "introduction"]], "Distributed-Memory Parallelization": [[8, "distributed-memory-parallelization"]], "Shared-Memory Parallelization": [[8, "shared-memory-parallelization"]], "Set Up": [[9, "set-up"]], "Anaconda Installation": [[9, "anaconda-installation"]], "Windows": [[9, "windows"]], "WSL": [[9, "wsl"]], "Without WSL": [[9, "without-wsl"]], "Mac OS and Linux": [[9, "mac-os-and-linux"]], "Creating a Parallel Programming conda environment": [[9, "creating-a-parallel-programming-conda-environment"]], "Confirming Install": [[9, "confirming-install"]]}, "indexentries": {}})