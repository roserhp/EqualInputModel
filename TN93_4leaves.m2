restart
K=frac(QQ[p_A,p_C,p_G,p_T]);
R=K[b_1,c_1,e_1,b_2,c_2,e_2,b_3,c_3,e_3,b_4,c_4,e_4,b_5,c_5,e_5];
i=1;
m=matrix{{0,b_i,c_i,b_i},{b_i,0,b_i,e_i},{c_i,b_i,0,b_i},{b_i,e_i,b_i,0}}
OE={A,C,G,T}
mm=mutableMatrix {toList apply(OE,i->p_(OE_((position(OE,j->j==i))%4))*(matrix m)_{position(OE,j->j==i)})}
for i to 3 do mm_(i,i)=1-sum flatten entries (matrix mm)^{i}
M1=matrix mm
i=2
M2=sub(M1,{b_1=>b_i,c_1=>c_i,e_1=>e_i})
i=3
M3=sub(M1,{b_1=>b_i,c_1=>c_i,e_1=>e_i})
i=4
M4=sub(M1,{b_1=>b_i,c_1=>c_i,e_1=>e_i})
i=5
M5=sub(M1,{b_1=>b_i,c_1=>c_i,e_1=>e_i})

--identidad en las hojas
M1=id_(R^4)
M2=M1
M3=M1
M4=M1

S=sort elements (set OE)^**4/splice/splice

pp=mutableMatrix(R,256,1)

for i to 255 do (
	pp_(i,0)=sum flatten toList apply(OE,k->apply(OE,kk->p_k*M1_(position(OE,l->l==k),position(OE,l->l==(S_i)_0))*M2_(position(OE,l->l==k),position(OE,l->l==(S_i)_1))*M3_(position(OE,l->l==kk),position(OE,l->l==(S_i)_2))*M4_(position(OE,l->l==kk),position(OE,l->l==(S_i)_3))*M5_(position(OE,l->l==k),position(OE,l->l==kk))))
) 

pp=matrix pp;

--ORIGINAL BASIS
--Repeated entries?
nonrepeatedpp=time unique(flatten entries pp);
length nonrepeatedpp --214 --> 256-214=42 equalities

pp_(position(S,i->i==(A,C,G,T)),0)==pp_(position(S,i->i==(A,T,G,C)),0)--true
pp_(position(S,i->i==(C,G,T,A)),0)==pp_(position(S,i->i==(T,G,C,A)),0)--true
pp_(position(S,i->i==(C,A,G,T)),0)==pp_(position(S,i->i==(T,A,G,C)),0)--true

pp_(position(S,i->i==(C,C,G,T)),0)==pp_(position(S,i->i==(T,C,G,C)),0)--false
pp_(position(S,i->i==(C,G,G,T)),0)==pp_(position(S,i->i==(T,G,G,C)),0)--true
pp_(position(S,i->i==(C,T,G,T)),0)==pp_(position(S,i->i==(T,T,G,C)),0)--false

p_G*pp_(position(S,i->i==(A,A,C,T)),0)==p_A*pp_(position(S,i->i==(A,G,C,T)),0)--false

equalEntries=select(unique for i to 255 list S_(positions(flatten entries pp,l->l==pp_(i,0))),k->length k>1);
length equalEntries
netList equalEntries
"4leavesOriginalBasisEqualEntries.txt" << netList equalEntries << endl << close

-- NEW BASIS
H=transpose(matrix{{4,4,4,4},{0,1/p_C,0,-1/p_T},{1/(p_A+p_G),-1/(p_C+p_T),1/(p_A+p_G),-1/(p_C+p_T)},{1/p_A,0,-1/p_G,0}});
H4=(transpose H)**(transpose H)**(transpose H)**(transpose H);
pbar=time H4*pp;

length select(flatten entries pbar,i->i==0) --176
zeroEntries=S_(positions(flatten entries pbar,i->i==0))

"4leavesBasisHzeroEntries.txt" << netList zeroEntries << endl << close

nonZeroEntries=S_(positions(flatten entries pbar,i->i!=0))
length unique nonZeroEntries --80

length select(zeroEntries,i->i_0==A and i_1==A) --12
length select(zeroEntries,i->i_0==A and i_1==C) --11
length select(zeroEntries,i->i_0==A and i_1==G) --11
length select(zeroEntries,i->i_0==A and i_1==T) --11
netList select(nonZeroEntries,i->i_0==A and i_1==T)


length select(zeroEntries,i->i_0==C and i_1==A) --11
length select(zeroEntries,i->i_0==C and i_1==C) --6
netList select(zeroEntries,i->i_0==C and i_1==C) 
length select(zeroEntries,i->i_0==C and i_1==G) --11
length select(zeroEntries,i->i_0==C and i_1==T) --16


length select(zeroEntries,i->i_0==G and i_1==A) --11
length select(zeroEntries,i->i_0==G and i_1==C) --11
length select(zeroEntries,i->i_0==G and i_1==G) --10
length select(zeroEntries,i->i_0==G and i_1==T) --11


length select(zeroEntries,i->i_0==T and i_1==A) --11
length select(zeroEntries,i->i_0==T and i_1==C) --16
length select(zeroEntries,i->i_0==T and i_1==G) --11
length select(zeroEntries,i->i_0==T and i_1==T) --6
netList select(zeroEntries,i->i_0==T and i_1==T) 


--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------


restart
K=frac(QQ[p_A,p_C,p_G,p_T]);
Rl=K[l_(5,A),l_(5,C),l_(5,G),l_(5,T)]
Hl=transpose(matrix{{1,1,1,1},{0,1/p_C,0,-1/p_T},{1/(p_A+p_G),-1/(p_C+p_T),1/(p_A+p_G),-1/(p_C+p_T)},{1/p_A,0,-1/p_G,0}});
--sanity check
Hl
inverse (transpose Hl)
Ml=Hl*diagonalMatrix(Rl,4,4,{l_(5,A),l_(5,C),l_(5,G),l_(5,T)})*inverse(Hl)
--sanity check
inverse(Hl)*Ml*Hl

--identidad en las hojas
M1=id_(Rl^4)
M2=M1
M3=M1
M4=M1

OE={A,C,G,T}

S=sort elements (set OE)^**4/splice/splice
netList S

qq=mutableMatrix(Rl,256,1)

for i to 255 do (
	qq_(i,0)=sum flatten toList apply(OE,k->apply(OE,kk->p_k*M1_(position(OE,l->l==k),position(OE,l->l==(S_i)_0))*M2_(position(OE,l->l==k),position(OE,l->l==(S_i)_1))*M3_(position(OE,l->l==kk),position(OE,l->l==(S_i)_2))*M4_(position(OE,l->l==kk),position(OE,l->l==(S_i)_3))*Ml_(position(OE,l->l==k),position(OE,l->l==kk))))
) 

qq=matrix qq;
netList (flatten entries qq)

H4=(transpose Hl)**(transpose Hl)**(transpose Hl)**(transpose Hl);

qbar=time H4*qq;
"4leaves_tensor_id_basisH.txt" << toString qbar << endl << close

netList toList apply(0..255,i->(S_i,qbar_(i,0)))
--only non-monomial entries
nonMonomial=select(S,i->(length terms qbar_(position(S,j->j==i),0)>1))
length nonMonomial
nonZeroEntries=S_(positions(flatten entries qbar,i->i!=0))
length nonZeroEntries 
monomialNonZeroEntries=select(nonZeroEntries,i->(length terms qbar_(position(S,j->j==i),0)==1))
length monomialNonZeroEntries

--Flattening 12|34 for identity at the leaves
flattql=mutableMatrix(Rl,16,16)
for i to 15 do (for j to 15 do 
    (if(j%4==j//4 and i%4==i//4) then flattql_(i,j)=p_(OE_(i%4))*Ml_(i%4,j%4));      
    );
flattql=matrix(flattql);
flattql=matrix(flattql);
--Change of basis
flattQl=time (transpose(Hl)**transpose(Hl))*flattql*(Hl**Hl);
--Quasi-block form
blockQl=flattQl_{0,3,12,14,11,15,10,2,8,5,1,4,9,6,13,7}^{0,3,12,14,11,15,10,2,8,5,1,4,9,6,13,7};
rank blockQl_{1,2,3,4}--1
rank blockQl_{14,15}--0
rank blockQl_{10,11,12,13}--1

--Flattening 12|34 general
Rgeneral=K[l_(1,A),l_(1,C),l_(1,G),l_(1,T),l_(2,A),l_(2,C),l_(2,G),l_(2,T),l_(3,A),l_(3,C),l_(3,G),l_(3,T),l_(4,A),l_(4,C),l_(4,G),l_(4,T),l_(5,A),l_(5,C),l_(5,G),l_(5,T)]
s={(A,A),(A,T),(T,A),(T,G),(G,T),(T,T),(G,G),(A,G),(G,A),(C,C),(A,C),(C,A),(C,G),(G,C),(T,C),(C,T)}
flattQl=sub(blockQl,Rgeneral);
flattPl=matrix toList apply(0..15,i->toList apply(0..15,j->l_(1,(s_i)_0)*l_(2,(s_i)_1)*l_(3,(s_j)_0)*l_(4,(s_j)_1)*flattQl_(i,j)));
--rk 0
rank flattPl_{position(s,i->i==(C,T)),position(s,i->i==(T,C))}--0
--rk 1
rank flattPl_{position(s,i->i==(A,T)),position(s,i->i==(T,A)),position(s,i->i==(G,T)),position(s,i->i==(T,G))}--1
rank flattPl_{position(s,i->i==(A,C)),position(s,i->i==(C,A)),position(s,i->i==(G,C)),position(s,i->i==(C,G))}--1
rank flattPl_{position(s,i->i==(A,G)),position(s,i->i==(G,A))}--1
--rk 2
rank flattPl_{position(s,i->i==(A,A)),position(s,i->i==(A,G)),position(s,i->i==(G,G))}--3
-- con variaciones
--rk 3
rank flattPl_{position(s,i->i==(A,A)),position(s,i->i==(A,G)),position(s,i->i==(A,T)),position(s,i->i==(T,T))}--3
rank flattPl_{position(s,i->i==(A,A)),position(s,i->i==(A,G)),position(s,i->i==(A,C)),position(s,i->i==(C,C))}--3
--con variaciones
--rk 4: all
rank flattPl_{position(s,i->i==(A,A)),position(s,i->i==(A,G)),position(s,i->i==(A,C)),position(s,i->i==(T,T))}--4

flattPl_{position(s,i->i==(A,A)),position(s,i->i==(A,C)),position(s,i->i==(A,G)),position(s,i->i==(A,T))}^{position(s,i->i==(A,A)),position(s,i->i==(A,C)),position(s,i->i==(A,G)),position(s,i->i==(A,T))}



syz flattPl --28
s_0,s_6,s_8
kernelPl=ker flattPl;
rank kernelPl
relations kernelPl --0
isFreeModule kernelPl --false
numgens kernelPl --28
mingens kernelPl --28
syz flattQl --12
kernelQl=ker flattQl;
rank kernelQl --12
relations kernelQl --0
isFreeModule kernelQl --false
numgens kernelQl --12
--------------
QBAR=sub(pbar^(positions(flatten entries pbar,i->i!=0)),Rgeneral);

PBAR=toList apply(nonZeroEntries,i->l_(1,i_0)*l_(2,i_1)*l_(3,i_2)*l_(4,i_3)*QBAR_(position(nonZeroEntries,j->j==i),0));

netList PBAR

varp=toList apply(nonZeroEntries,i->(symbol p)_i);
S=K[varp]; 
f=map(Rgeneral,S,PBAR);
f(p_(A,A,A,A))
f(p_(A,T,T,T))
I=time kernel f;
betti I
J=trim I;
betti J

--Try to compute first the monomial parametrization: get rid of combinations of CC,GG,TT
toricPBAR=toList apply(monomialNonZeroEntries,i->l_(1,i_0)*l_(2,i_1)*l_(3,i_2)*l_(4,i_3)*QBAR_(position(nonZeroEntries,j->j==i),0));
netList toricPBAR


toricvarp=toList apply(monomialNonZeroEntries,i->(symbol p)_i);
toricS=K[toricvarp]; 
toricf=map(Rgeneral,toricS,toricPBAR);
toricf(p_(A,A,A,A))
toricf(p_(A,T,T,T))
toricI=time kernel toricf;
Too many heap sections: Increase MAXHINCR or MAX_HEAP_SECTS
Aborted (core dumped)

-------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------

--------------------------------------------------------------------
--------------------------------------------------------------------
-------------------------------------------------------------------

restart
K=frac(QQ[p_A,p_C,p_G,p_T]);
Rl=K[l_(5,A),l_(5,C),l_(5,G),l_(5,T)]
qbar=value get "4leaves_tensor_id_basisH.txt";

OE={A,C,G,T}
S=sort elements (set OE)^**4/splice/splice;

--Flattening 12|34
s={(A,A),(A,T),(T,A),(T,G),(G,T),(T,T),(G,G),(A,G),(G,A),(C,C),(A,C),(C,A),(C,G),(G,C),(T,C),(C,T)}
flattq=mutableMatrix(Rl,4^2,4^2);
for i to 15 do (for j to 15 do
    flattq_(i,j)=qbar_(position(S,k->k==((s_i)_0,(s_i)_1,(s_j)_0,(s_j)_1)),0);
    );
flattq=matrix flattq;

Rgeneral=K[l_(1,A),l_(1,C),l_(1,G),l_(1,T),l_(2,A),l_(2,C),l_(2,G),l_(2,T),l_(3,A),l_(3,C),l_(3,G),l_(3,T),l_(4,A),l_(4,C),l_(4,G),l_(4,T),l_(5,A),l_(5,C),l_(5,G),l_(5,T)]
qbar=sub(qbar,Rgeneral);
pbar=transpose matrix{toList apply(S,i->l_(1,i_0)*l_(2,i_1)*l_(3,i_2)*l_(4,i_3)*qbar_(position(S,j->j==i),0))};
flattq=sub(flattq,Rgeneral);
flattp=matrix toList apply(0..15,i->toList apply(0..15,j->l_(1,(s_i)_0)*l_(2,(s_i)_1)*l_(3,(s_j)_0)*l_(4,(s_j)_1)*flattq_(i,j)));

flattp_{position(s,i->i==(A,A)),position(s,i->i==(A,C)),position(s,i->i==(A,G)),position(s,i->i==(A,T))}^{position(s,i->i==(A,A)),position(s,i->i==(A,C)),position(s,i->i==(A,G)),position(s,i->i==(A,T))}

var=toList apply(S,i->(symbol x)_i)
Rx=K[var]
f=map(Rgeneral,Rx,flatten entries pbar);
--Example
f(x_(A,A,A,A))
f(x_(A,A,A,G))

R3=K[toList apply({(A,A,A),(A,C,C),(A,G,G),(A,T,T),(C,A,C),(C,C,A),(C,C,C),(C,C,G),(C,G,C),(G,A,G),(G,C,C),(G,G,A),(G,G,G),(G,T,T),(T,A,T),(T,G,T),(T,T,A),(T,T,G),(T,T,T)},i->(symbol p)_i)]
I=value get "3leavesVanishingIdeal.txt";
CI=ideal{I_0,I_1,I_2,I_3,I_4,I_6,I_13,I_28,I_37};
netList CI_*

s4=apply({(A,A,A),(A,C,C),(A,G,G),(A,T,T),(C,A,C),(C,C,A),(C,C,C),(C,C,G),(C,G,C),(G,A,G),(G,C,C),(G,G,A),(G,G,G),(G,T,T),(T,A,T),(T,G,T),(T,T,A),(T,T,G),(T,T,T)},i->append(i,A));
s1=apply({(A,A,A),(A,C,C),(A,G,G),(A,T,T),(C,A,C),(C,C,A),(C,C,C),(C,C,G),(C,G,C),(G,A,G),(G,C,C),(G,G,A),(G,G,G),(G,T,T),(T,A,T),(T,G,T),(T,T,A),(T,T,G),(T,T,T)},i-> sequence A|i);

g4=map(Rx,R3,toList apply(s4,i->x_i));
CI4=g4(CI);
netList CI4_*
apply(flatten entries gens CI4,i->f(i))

g1=map(Rx,R3,toList apply(s1,i->x_i));
CI1=g1(CI);
netList CI1_*
apply(flatten entries gens CI1,i->f(i))

CItripod=CI1+CI4; --is this really a model invariant??? 
--it has to be because adding A doesn't change anything
--if the split was 14|23 then those two would be in the same side...maybe it wouldn't be giving a local complete intersection then?
time codim CItripod --14 => not a complete intersection
-- used 234.656 seconds
time codim CI1 --9
time codim CI4 --9
betti CItripod --18
CItripod==trim CItripod --true

s={(A,A),(A,T),(T,A),(T,G),(G,T),(T,T),(G,G),(A,G),(G,A),(C,C),(A,C),(C,A),(C,G),(G,C),(T,C),(C,T)}

flattp_{position(s,i->i==(A,T)),position(s,i->i==(T,A)),position(s,i->i==(T,G)),position(s,i->i==(G,T))}
rank oo
R1T=flattp_{position(s,i->i==(A,T)),position(s,i->i==(T,A)),position(s,i->i==(T,G)),position(s,i->i==(G,T))}^{position(s,i->i==(A,T)),position(s,i->i==(T,A)),position(s,i->i==(T,G)),position(s,i->i==(G,T)),position(s,i->i==(T,T))}

         AT         j
AT    (A,T,A,T)   (AT,j)
i        (i,AT)    (i,j)

rk1T=flatten toList apply({(T,A),(T,G),(G,T)},j->apply({(T,A),(T,G),(G,T),(T,T)},i->x_(A,T,A,T)*x_(i|j)-x_((A,T)|j)*x_(i|(A,T))))
length rk1T --12

R1C=flattp_{position(s,i->i==(A,C)),position(s,i->i==(C,A)),position(s,i->i==(C,G)),position(s,i->i==(G,C))}^{position(s,i->i==(C,C)),position(s,i->i==(A,C)),position(s,i->i==(C,A)),position(s,i->i==(G,C)),position(s,i->i==(C,G))}

         AC         j
AC    (A,C,A,C)   (AC,j)
i        (i,AC)    (i,j)

rk1C=flatten toList apply({(C,A),(C,G),(G,C)},j->apply({(C,A),(C,G),(G,C),(C,C)},i->x_(A,C,A,C)*x_(i|j)-x_((A,C)|j)*x_(i|(A,C))))
length rk1C --12

apply(rk1C,i->f(i))
apply(rk1T,i->f(i))
f(x_(A,C,C,A)*x_(C,A,A,C))
f(x_(A,C,A,C)*x_(C,A,C,A))


flattp_{position(s,i->i==(A,G)),position(s,i->i==(G,A))}
rank oo
flattp_{position(s,i->i==(A,G)),position(s,i->i==(G,A))}^{position(s,i->i==(T,T)),position(s,i->i==(G,G)),position(s,i->i==(A,G)),position(s,i->i==(G,A)),position(s,i->i==(C,C))}
         AG         GA
AG    (A,G,A,G)   (A,G,G,A)
i        (i,AG)    (i,GA)
rk1G=toList apply({(G,A),(C,C),(G,G),(T,T)},i->x_(A,G,A,G)*x_(i|(G,A))-x_(A,G,G,A)*x_(i|(A,G)))
length rk1G --4
apply(rk1G,i->f(i))

rk1=rk1C|rk1G|rk1T
length rk1

18
12*2+4

R2G=flattp_{position(s,i->i==(A,A)),position(s,i->i==(A,G)),position(s,i->i==(G,G))}^{position(s,i->i==(A,A)),position(s,i->i==(T,T)),position(s,i->i==(G,G)),position(s,i->i==(A,G)),position(s,i->i==(G,A)),position(s,i->i==(C,C))}

flattp_{position(s,i->i==(A,A)),position(s,i->i==(A,G))}^{position(s,i->i==(A,A)),position(s,i->i==(A,G))}


        AA     AG     GG
   
AA     AAAA   AAAG   AAGG

AG     AGAA   AGAG   AGGG

i      iAA     iAG   iGG


rk2=toList apply({(G,A),(C,C),(G,G),(T,T)},i->det matrix{{x_(A,A,A,A),x_(A,A,A,G),x_(A,A,G,G)},
	{x_(A,G,A,A),x_(A,G,A,G),x_(A,G,G,G)},{x_(i|(A,A)),x_(i|(A,G)),x_(i|(G,G))}})

apply(rk2,i->f(i))

length rk2 --4

--rk 3 T
flattp_{position(s,i->i==(A,A)),position(s,i->i==(A,G)),position(s,i->i==(A,T)),position(s,i->i==(T,T))}
rank oo

0,1,2,3,4,5,6,7,8,9
length{(A,A),(A,T),(T,A),(T,G),(G,T),(T,T),(G,G),(A,G),(G,A),(C,C)}


        AA     AG     AT   TT
   
AA     AAAA   AAAG   AAAT  AATT

AG     AGAA   AGAG   AGAT  AGTT

AT     ATAA   ATAG   ATAT  ATTT

i      iAA     iAG   iAT   iTT

rk3T=toList apply({(T,A),(T,G),(G,T),(T,T),(G,G),(G,A),(C,C)},i->det matrix{{x_(A,A,A,A),x_(A,A,A,G),x_(A,A,A,T),x_(A,A,T,T)},
	{x_(A,G,A,A),x_(A,G,A,G),x_(A,G,A,T),x_(A,G,T,T)},{x_(A,T,A,A),x_(A,T,A,G),x_(A,T,A,T),x_(A,T,T,T)},
	{x_(i|(A,A)),x_(i|(A,G)),x_(i|(A,T)),x_(i|(T,T))}})

apply(rk3T,i->f(i))
length rk3T --7


--rk 3 C
flattp_{position(s,i->i==(A,A)),position(s,i->i==(A,G)),position(s,i->i==(A,C)),position(s,i->i==(C,C))}

0,5,6,7,8,9,10,11,12,13
length{(A,A),(T,T),(G,G),(A,G),(G,A),(C,C),(A,C),(C,A),(C,G),(G,C)}


        AA     AG     AC   CC
   
AA     AAAA   AAAG   AAAC  AACC

AG     AGAA   AGAG   AGAC  AGCC

AC     ACAA   ACAG   ACAC  ACCC

i      iAA     iAG   iAC   iCC

rk3C=toList apply({(T,T),(G,G),(G,A),(C,C),(C,A),(C,G),(G,C)},i->
    det matrix{{x_(A,A,A,A),x_(A,A,A,G),x_(A,A,A,C),x_(A,A,C,C)},
	{x_(A,G,A,A),x_(A,G,A,G),x_(A,G,A,C),x_(A,G,C,C)},{x_(A,C,A,A),x_(A,C,A,G),x_(A,C,A,C),x_(A,C,C,C)},
	{x_(i|(A,A)),x_(i|(A,G)),x_(i|(A,C)),x_(i|(C,C))}})

apply(rk3C,i->f(i))
length rk3C


rk3=rk3C|rk3T
length rk3

rk0={x_(C,T,C,T),x_(C,T,T,C),x_(T,C,C,T),x_(T,C,T,C)}
apply(rk0,i->f(i))

edge=rk0|rk1|rk2|rk3
CIedge=ideal edge;
CI=CIedge+CItripod;
betti CI --68
numgens CI  --18 
numgens CIedge --50

jac=jacobian CI;
time rank jac --68 What is this exactly computing?

ident=flatten entries sub(sub(sub(qbar,Rl),apply(gens Rl,i->i=>1)),Rx);
jacId=sub(jac,matrix{ident});
time rank jacId --68
