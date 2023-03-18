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

--Flattening 12|34 for identity at the leaves
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

restart
K=frac(QQ[p_A,p_C,p_G,p_T]);
H=transpose(matrix{{1,1,1,1},{0,1/p_C,0,-1/p_T},{1/(p_A+p_G),-1/(p_C+p_T),1/(p_A+p_G),-1/(p_C+p_T)},{1/p_A,0,-1/p_G,0}});
Hpi=sub(H,{sub(p_A,K)=>1/2,sub(p_C,K)=>1/3,sub(p_G,K)=>1/8,sub(p_T,K)=>1/24})
Ht1=inverse transpose Hpi
D=diagonalMatrix({1/2,1/3,1/8,1/24})
G=inverse D
for i to 3 do for j from i to 3 do print (i,j,(transpose(Ht1_{i})*G*Ht1_{j})_(0,0))
