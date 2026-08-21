import re

def parse_ln(ln):
	ln_list = []
	for s in ln:
		if s != "\n":
			ln_list.append(F(s))
	return ln_list

q=7
b=2
H=[]
L=[]
T=[]
MW_Calls=[]
n_list=[32,64,128,256,512,1024]
for n in n_list:
	H_row = []
	L_row = []
	T_row = []
	MW_Calls_row = []
	# Algorithm 0
	stat_file = "experiments_q" + str(q) + "/LeeBKZ_q" + str(q) + "_b" + str(b) + "_n" + str(n) + "_0.out" # Hamming BKZ
	f = open(stat_file, "r")
	ln=''
	while 'Reduced basis' not in ln:
		ln=f.readline()
	red_basis = re.findall(r"[-+]?(?:\d*\.*\d+)", ln)
	H_row.append(red_basis[1])
	L_row.append(red_basis[3])
	ln=f.readline() # Avg.time
	time = re.findall(r"[-+]?(?:\d*\.*\d+)", ln)
	T_row.append(time[0])
	ln=f.readline() # Avg.Mw_codeword.calls
	mwc = re.findall(r"[-+]?(?:\d*\.*\d+)", ln)
	MW_Calls_row.append(mwc[0])
	f.close()

	# Algorithm 1
	stat_file = "experiments_q" + str(q) + "/LeeBKZ_q" + str(q) + "_b" + str(b) + "_n" + str(n) + "_1.out" # Hamming BKZ selecting lowest Lee wt.
	f = open(stat_file, "r")
	ln=''
	while 'Reduced basis' not in ln:
		ln=f.readline()
	red_basis = re.findall(r"[-+]?(?:\d*\.*\d+)", ln)
	H_row.append(red_basis[1])
	L_row.append(red_basis[3])
	ln=f.readline() # Avg.time
	time = re.findall(r"[-+]?(?:\d*\.*\d+)", ln)
	T_row.append(time[0])
	ln=f.readline() # Avg.Mw_codeword.calls
	mwc = re.findall(r"[-+]?(?:\d*\.*\d+)", ln)
	MW_Calls_row.append(mwc[0])
	f.close()

	# Algorithm 2
	stat_file = "experiments_q" + str(q) + "/LeeBKZ_q" + str(q) + "_b" + str(b) + "_n" + str(n) + "_2.out" # BKZ with minimum LeeWt. codeword + primitive
	f = open(stat_file, "r")
	ln=''
	while 'Reduced basis' not in ln:
		ln=f.readline()
	red_basis = re.findall(r"[-+]?(?:\d*\.*\d+)", ln)
	H_row.append(red_basis[1])
	L_row.append(red_basis[3])
	ln=f.readline() # Avg.time
	time = re.findall(r"[-+]?(?:\d*\.*\d+)", ln)
	T_row.append(time[0])
	ln=f.readline() # Avg.Mw_codeword.calls
	mwc = re.findall(r"[-+]?(?:\d*\.*\d+)", ln)
	MW_Calls_row.append(mwc[0])
	f.close()

	# Algorithm 3 (LeeBKZ)
	stat_file = "experiments_q" + str(q) + "/LeeBKZ_q" + str(q) + "_b" + str(b) + "_n" + str(n) + "_3.out" # LeeBKZ
	f = open(stat_file, "r")
	ln=''
	while 'Reduced basis' not in ln:
		ln=f.readline()
	red_basis = re.findall(r"[-+]?(?:\d*\.*\d+)", ln)
	H_row.append(red_basis[1])
	L_row.append(red_basis[3])
	ln=f.readline() # Avg.time
	time = re.findall(r"[-+]?(?:\d*\.*\d+)", ln)
	T_row.append(time[0])
	ln=f.readline() # Avg.Mw_codeword.calls
	mwc = re.findall(r"[-+]?(?:\d*\.*\d+)", ln)
	MW_Calls_row.append(mwc[0])
	f.close()

	# Algorithm 3 (LLL)
	if b==2:
		stat_file = "experiments_q" + str(q) + "/LeeBKZ_q" + str(q) + "_b" + str(b) + "_n" + str(n) + "_LLL.out"
		f = open(stat_file, "r")
		ln=''
		while 'Reduced basis' not in ln:
			ln=f.readline()
		red_basis = re.findall(r"[-+]?(?:\d*\.*\d+)", ln)
		H_row.append(red_basis[1])
		L_row.append(red_basis[3])
		ln=f.readline() # Avg.time
		time = re.findall(r"[-+]?(?:\d*\.*\d+)", ln)
		T_row.append(time[0])
		ln=f.readline() # Avg.Mw_codeword.calls
		mwc = re.findall(r"[-+]?(?:\d*\.*\d+)", ln)
		MW_Calls_row.append(mwc[0])
		f.close()

	H.append(H_row)
	L.append(L_row)
	T.append(T_row)
	MW_Calls.append(MW_Calls_row)

# Output stats to terminal
print(f"q={q}, b={b}")
print("n | Alg0 | Alg1 | Alg2 | Alg3 (LeeBKZ) | Alg3 (LeeLLL, only for beta=2)")
print("Lee reduction quality (ell_1 in LEE Profile):")
i=0
for n in n_list:
	print(f"n={n} ", end="")
	print(*L[i])
	i=i+1

print("Hamming reduction quality (ell_1 in Hamming Profile):")
i=0
for n in n_list:
	print(f"n={n} ", end="")
	print(*H[i])
	i=i+1

print("Number of MinimumWeightCodeword() calls:")
i=0
for n in n_list:
	print(f"n={n} ", end="")
	print(*MW_Calls[i])
	i=i+1

print("Operation time, s.:")
i=0
for n in n_list:
	print(f"n={n} ", end="")
	print(*T[i])
	i=i+1
