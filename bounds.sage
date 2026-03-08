load("LEE.sage")

def LeeSphereSize(n,t,q):
    R.<z> = PolynomialRing(ZZ)
    Sq1 = 1
    for i in [1..floor((q-1)/2)]:
        Sq1 += 2*z^i
    if q%2 == 0:
        Sq1 += z^(q/2)
    return (Sq1^n).coefficients()[t]

def LeeBallSize(n,t,q):
    return sum(LeeSphereSize(n,r,q) for r in [0..t])
    
def GV_d_min_LEE(q,n,k):
    d=1
    while LeeBallSize(n,d-1,q)<q^(n-k):
        d=d+1
    return d

def GV_d_min_LEE1(q,n,k):
    R.<z> = PolynomialRing(ZZ)
    Sq1 = 1
    for i in [1..floor((q-1)/2)]:
        Sq1 += 2*z^i
    if q%2 == 0:
        Sq1 += z^(q/2)
    Sq1 = Sq1^n
    d=1
    LLeeBallSize=Sq1.coefficients()[d-1]
    while LLeeBallSize<q^(n-k):
        d=d+1
        LLeeBallSize = LLeeBallSize*2 + Sq1.coefficients()[d-1]
    return d

def GV(q,n,k):
    d=1;
    while sum(binomial(n,j)*(q-1)^j for j in [0..d-1]) < q^(n-k):
        d+=1;
    return d

for q in [5,11,19,79,127]:
    for n in [32,64,128,256,512,1024]:
        k=int(n/2)
        print(f"q={q} n={n} k={k} GV d_min LEE={GV_d_min_LEE1(q,n,k)}")
        print(f"q={q} n={n} k={k} GV d_min HAM={GV(q,n,k)}")
        print("---")
