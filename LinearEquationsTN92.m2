--TESTS

restart
K=frac(QQ[p_A,p_C,p_G,p_T])
R=K[l_(1,A),l_(1,C),l_(1,G),l_(1,T),l_(2,A),l_(2,C),l_(2,G),l_(2,T),l_(3,A),l_(3,C),l_(3,G),l_(3,T),l_(4,A),l_(4,C),l_(4,G),l_(4,T),l_(5,A),l_(5,C),l_(5,G),l_(5,T)]

--aaaa=256*l_(5,A)
atat=16*(1/p_A+1/p_G)*l_(5,T)
atgt=4/(p_A+p_G)*(1/p_A+1/p_G)*l_(5,T)
gtgt=1/(p_A+p_G)^2*(1/p_A+1/p_G)*l_(5,T)

qI=matrix{{atat,atat,atgt,atgt},{atat,atat,atgt,atgt},{atgt,atgt,gtgt,gtgt},{atgt,atgt,gtgt,gtgt}}

s={(A,T),(T,A),(T,G),(G,T)}

q=matrix toList apply(0..3,i->toList apply(0..3,j->l_(1,(s_i)_0)*l_(2,(s_i)_1)*l_(3,(s_j)_0)*l_(4,(s_j)_1)*qI_(i,j)))

varp=flatten toList apply(0..3,i->toList apply(0..3,j->(symbol p)_(s_i,s_j)))

S=QQ[varp]

P=transpose genericMatrix(S,p_((A,T),(A,T)),4,4)

RS=K[gens S,gens R]
P=sub(P,RS)
q=sub(q,RS)
I=ideal flatten entries (P-q);
aux=ideal toList apply(1..5,i->(p_A+p_G)*(l_(i,A)-l_(i,C))-(1-(p_A+p_G))*(l_(i,A)-l_(i,T))-(2*(p_A+p_G)-1)*(l_(i,A)-l_(i,G)))
Iaux=trim(I+aux);    
netList Iaux_*
Iaux==I

J=time eliminate(apply(gens R,i->sub(i,RS)),I);
betti (trim J)
netList J_*
Jaux=time eliminate(apply(gens R,i->sub(i,RS)),Iaux);
betti Jaux
netList Jaux_*
--Ring with lambda's first and eliminate 20


--Checks (correct parametrization in terms of eigenvalues?)
restart
K=frac(QQ[p_A,p_C,p_G,p_T])
R=K[b,c]
m=matrix{{0,b,c,b},{b,0,b,c},{c,b,0,b},{b,c,b,0}}
mm=mutableMatrix toList apply(0..3,i->(gens K)_i*(matrix m)_i)
for i to 3 do mm_(i,i)=1-sum flatten entries (matrix mm)^{i}
M=matrix mm
--Remark 4
1/(1/p_A+1/p_G)*((1/p_A)*M_(0,0)-(1/p_G)*M_(0,2)-(1/p_A)*M_(2,0)+(1/p_G)*M_(2,2))
--Check computation of q_ATTT
((p_A^2*p_G^2)/(p_G^2-p_A^2))*((1/p_A^2)*M_(0,0)+(1/p_G^2)*M_(0,2)-(1/p_A^2)*M_(2,0)-(1/p_G^2)*M_(2,2))
M_(2,0)+M_(2,2)-M_(0,0)-M_(0,2)
--Check computation of q_GCCC
((p_C^2*p_T^2)/(p_T^2-p_C^2))*((1/p_C^2)*M_(1,1)+(1/p_T^2)*M_(1,3)-(1/p_C^2)*M_(3,1)-(1/p_T^2)*M_(3,3))
M_(2,0)+M_(2,2)-M_(0,0)-M_(0,2)
--Check computation of q_GGTT
OE={A,C,G,T}
GGTT=sum apply(OE,i->(p_i/(p_i+p_(OE_((position(OE,j->j==i)+2)%4)))^2)*((1/p_A^2)*M_(position(OE,j->j==i),0)+(1/p_G^2)*M_(position(OE,j->j==i),2)))
netList terms GGTT
(((p_A+p_G)^2-(p_C+p_T)^2)/(p_A*p_G*(p_A+p_G)*(p_C+p_T)))*b==(terms GGTT)_0

-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------
-- ACTUAL COMPUTATIONS
-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------

restart
K=frac(QQ[p_A,p_C,p_G,p_T])
R=K[l_(1,A),l_(1,C),l_(1,G),l_(1,T),l_(2,A),l_(2,C),l_(2,G),l_(2,T),l_(3,A),l_(3,C),l_(3,G),l_(3,T),l_(4,A),l_(4,C),l_(4,G),l_(4,T),l_(5,A),l_(5,C),l_(5,G),l_(5,T)]

aaaa=256*l_(5,A)
aatt=16*(1/p_A+1/p_G)*l_(5,A)
aagg=16*(1/((p_A+p_G)*(p_C+p_T)))*l_(5,A)
aacc=16*(1/p_C+1/p_T)*l_(5,A)
atat=16*(1/p_A+1/p_G)*l_(5,T)
atgt=4/(p_A+p_G)*(1/p_A+1/p_G)*l_(5,T)
gtgt=1/(p_A+p_G)^2*(1/p_A+1/p_G)*l_(5,T)
attt=4*(p_G^2-p_A^2)/(p_A^2*p_G^2)*l_(5,T)
accc=4*(p_T^2-p_C^2)/(p_C^2*p_T^2)*l_(5,C)
gttt=(p_G-p_A)/(p_A^2*p_G^2)*l_(5,T)
gccc=(p_C-p_T)/(p_C^2*p_T^2)*l_(5,C)
tttt=(p_A+p_G)/(p_A^2*p_G^2)*(l_(5,A)-(p_C+p_T)*(l_(5,A)-l_(5,G))+(p_A-p_G)^2/(p_A*p_G)*l_(5,T))
ggtt=(1)/(p_A*p_G*(p_A+p_G)*(p_C+p_T))*((p_C+p_T-p_A-p_G)*l_(5,G)+(p_A+p_G)*l_(5,A))
ccgg=(1)/(p_C*p_T*(p_A+p_G)*(p_C+p_T))*((p_A+p_G-p_C-p_T)*l_(5,G)+(p_C+p_T)*l_(5,A))
cctt=(1/p_A+1/p_G)*(1/p_C+1/p_T)*(l_(5,A)-l_(5,G))
agtt=(4/(p_A*p_G))*l_(5,G)

gggg=(1-2*(p_C+p_T))^2/((p_A+p_G)^3*(p_C+p_T)^3)*(l_(5,G)-l_(5,A))+(1/(p_A+p_G)^3)*l_(5,A)+(1/(p_C+p_T)^3)*l_(5,A)

aggg=4*l_(5,G)*(1/(p_A+p_G)^2+1/(p_C+p_T)^2)
agag=16*l_(5,G)*(1/(p_A+p_G)+1/(p_C+p_T))
agcc=(-4/(p_C*p_T))*l_(5,G)

cccc=(p_C+p_T)/(p_C^2*p_T^2)*(l_(5,A)-(p_A+p_G)*(l_(5,A)-l_(5,G))+(p_C-p_T)^2/(p_C*p_T)*l_(5,C))


accc=4*(p_T^2-p_C^2)/(p_C^2*p_T^2)*l_(5,C)
gccc=(p_C-p_T)/(p_C^2*p_T^2)*l_(5,C)
acac=16*(1/p_C+1/p_T)*l_(5,C)
accg=-(4/(p_C+p_T))*(1/p_C+1/p_T)*l_(5,C)
cgcg=(1/(p_C+p_T)^2)*(1/p_C+1/p_T)*l_(5,C)

qI=matrix{
    {aaaa,0,0,0,0,aatt,aagg,0,0,aacc,0,0,0,0},
    {0,atat,atat,atgt,atgt,attt,0,0,0,0,0,0,0,0},
    {0,atat,atat,atgt,atgt,attt,0,0,0,0,0,0,0,0},
    {0,atgt,atgt,gtgt,gtgt,gttt,0,0,0,0,0,0,0,0},
    {0,atgt,atgt,gtgt,gtgt,gttt,0,0,0,0,0,0,0,0},
    {aatt,attt,attt,gttt,gttt,tttt,ggtt,agtt,agtt,cctt,0,0,0,0},
    {aagg,0,0,0,0,ggtt,gggg,aggg,aggg,ccgg,0,0,0,0},
       {0,0,0,0,0,agtt,aggg,agag,agag,agcc,0,0,0,0},
       {0,0,0,0,0,agtt,aggg,agag,agag,agcc,0,0,0,0},
    {aacc,0,0,0,0,cctt,ccgg,agcc,agcc,cccc,accc,accc,gccc,gccc},
                   {0,0,0,0,0,0,0,0,0,accc,acac,acac,accg,accg},
                   {0,0,0,0,0,0,0,0,0,accc,acac,acac,accg,accg},
                   {0,0,0,0,0,0,0,0,0,gccc,accg,accg,cgcg,cgcg},
                   {0,0,0,0,0,0,0,0,0,gccc,accg,accg,cgcg,cgcg}}



s={(A,A),(A,T),(T,A),(T,G),(G,T),(T,T),(G,G),(A,G),(G,A),(C,C),(A,C),(C,A),(C,G),(G,C)}
length s

q=matrix toList apply(0..13,i->toList apply(0..13,j->l_(1,(s_i)_0)*l_(2,(s_i)_1)*l_(3,(s_j)_0)*l_(4,(s_j)_1)*qI_(i,j)))

varp=flatten toList apply(0..13,i->toList apply(0..13,j->(symbol p)_(s_i,s_j)))

S=QQ[varp]

P=transpose genericMatrix(S,p_((A,A),(A,A)),14,14)

RS=K[gens S,gens R]
P=sub(P,RS)
q=sub(q,RS)
I=trim ideal flatten entries (P-q);
betti I
aux=ideal toList apply(1..5,i->(p_A+p_G)*(l_(i,A)-l_(i,C))-(1-(p_A+p_G))*(l_(i,A)-l_(i,T))-(2*(p_A+p_G)-1)*(l_(i,A)-l_(i,G)))
netList aux_*
-- 116 linear equations, 80 quintics
i=1
(p_A+p_G)*(l_(i,A)-l_(i,C))-(1-(p_A+p_G))*(l_(i,A)-l_(i,T))-(2*(p_A+p_G)-1)*(l_(i,A)-l_(i,G))
--l_(i,A) disappears from the equation
Iaux=trim(I+aux);    
Iaux==I
betti Iaux

J=time eliminate(apply(gens R,i->sub(i,RS)),I);
betti (trim J)
netList J_*
Jaux=time eliminate(apply(gens R,i->sub(i,RS)),Iaux);
betti Jaux
netList Jaux_*
