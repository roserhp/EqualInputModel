restart
K=frac(QQ[p_A,p_C,p_G,p_T]);
R=K[b_1,c_1,e_1,b_2,c_2,e_2,b_3,c_3,e_3];
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

--Is the product of two TN93 matrices a TN93 matrix with different parameters?
M=M1*M2

M_(1,0)
M1_(1,0)
1/p_T*M_(1,3)
M1_(1,3)
M1_(3,1)
1/p_C*M_(3,1)
1/p_C*M_(3,1)==1/p_T*M_(1,3)

--Identity on the leaves
restart
K=frac(QQ[p_A,p_C,p_G,p_T]);
R=K[b_1,c_1,e_1,b_2,c_2,e_2,b_3,c_3,e_3];
OE={A,C,G,T}
M1=id_(R^4)
M2=M1
M3=M1

S=sort elements (set OE)^**3/splice
netList S

pp=mutableMatrix(R,64,1)
for i to 63 do (
	pp_(i,0)=sum apply(OE,k->p_k*M1_(position(OE,l->l==k),position(OE,l->l==(S_i)_0))
	                           *M2_(position(OE,l->l==k),position(OE,l->l==(S_i)_1))
		                   *M3_(position(OE,l->l==k),position(OE,l->l==(S_i)_2)))
) 

pp=matrix pp;

--ORIGINAL BASIS
length select(flatten entries pp,i->i==0) --0
length unique(flatten entries pp) --52 --> 64-52=12 equalities

equalEntries=select(unique for i to 63 list S_(positions(flatten entries pp,l->l==pp_(i,0))),k->length k>1);
length equalEntries --12
netList equalEntries
"3leavesOriginalBasisEqualEntries.txt" << netList equalEntries << endl << close


--NEW BASIS
H=transpose(matrix{{4,4,4,4},{0,1/p_C,0,-1/p_T},{1/(p_A+p_G),-1/(p_C+p_T),1/(p_A+p_G),-1/(p_C+p_T)},{1/p_A,0,-1/p_G,0}});

H3=(transpose H)**(transpose H)**(transpose H);

pbar=time H3*pp;

nonZeroEntries=S_(positions(flatten entries pbar,i->i!=0))
length unique nonZeroEntries --19
zeroEntries=S_(positions(flatten entries pbar,i->i==0))
length zeroEntries --45

"3leavesBasisHzeroEntries.txt" << netList zeroEntries << endl << close

Rgeneral=K[l_(1,A),l_(1,C),l_(1,G),l_(1,T),l_(2,A),l_(2,C),l_(2,G),l_(2,T),l_(3,A),l_(3,C),l_(3,G),l_(3,T)]
QBAR=sub(pbar^(positions(flatten entries pbar,i->i!=0)),Rgeneral);
PBAR=toList apply(nonZeroEntries,i->l_(1,i_0)*l_(2,i_1)*l_(3,i_2)*QBAR_(position(nonZeroEntries,j->j==i),0));

netList PBAR

varp=toList apply(nonZeroEntries,i->(symbol p)_i);
S=K[varp]; --QQ or K??
f=map(Rgeneral,S,PBAR);
f(p_(A,A,A))
f(p_(T,T,T))
test=time trim kernel f;
--16 sec aprox
betti test
netList test_*
codim test
dim test
I2=ideal toList apply(0..8,i->test_i);
netList I2_*
codim I2
CI=time regSeqInIdeal(test,9);


apply(flatten entries gens test,i->length terms i)
------------------------------------------------------------------------------
-------------------------------------------------------------------------------
--Checks with original basis
--uniform distribution
uni=sub(pp,{p_A=>1/4,p_C=>1/4,p_G=>1/4,p_T=>1/4});
--b=c=e
EI=sub(pp,flatten toList apply(1..3,i->{c_i=>b_i,e_i=>b_i}));
nonrepeatedEI=time unique(flatten entries EI);
length nonrepeatedEI --44=64-20 
--matches Example pg 1119 in Casanellas-Steel
--but we need 16 more equations!!!
p_T*EI_(position(S,i->i==(A,A,G)),0)==p_G*EI_(position(S,i->i==(A,A,T)),0) 
p_T*EI_(position(S,i->i==(A,A,A)),0)==p_G*EI_(position(S,i->i==(A,A,T)),0) --false
p_T*EI_(position(S,i->i==(A,G,A)),0)==p_G*EI_(position(S,i->i==(A,T,A)),0) 
p_T*EI_(position(S,i->i==(G,A,A)),0)==p_G*EI_(position(S,i->i==(T,A,A)),0) 
p_T*EI_(position(S,i->i==(A,C,G)),0)==p_G*EI_(position(S,i->i==(A,C,T)),0) 
EI_(position(S,i->i==(A,C,G)),0)==EI_(position(S,i->i==(G,A,C)),0) 
p_C*p_T*EI_(position(S,i->i==(A,A,G)),0)==p_A*p_A*EI_(position(S,i->i==(C,G,T)),0)--false


--b=c=e and uniform distr
JC=sub(uni,flatten toList apply(1..3,i->{c_i=>b_i,e_i=>b_i}));
nonrepeatedJC=time unique(flatten entries JC);
length nonrepeatedJC --5=64-59 matches Casanellas-Fernandez-Kedzierska
JC_(position(S,i->i==(A,A,A)),0)==JC_(position(S,i->i==(T,T,T)),0) --true
JC_(position(S,i->i==(A,A,G)),0)==JC_(position(S,i->i==(A,A,T)),0) --true

--c=e and uniform distr
K80=sub(uni,flatten toList apply(1..3,i->{e_i=>c_i}));
nonrepeatedK80=time unique(flatten entries K80);
length nonrepeatedK80 --10=64-54 matches Casanellas-Fernandez-Kedzierska

----------------------------------------------------------------------------
----------------------------------------------------------------------------
-- Checks with original basis rescaled with inverse of pi's appearing in the configuration
qq=mutableMatrix pp;
for i to 63 do (
	qq_(i,0)=1/(p_((S_i)_0)*p_((S_i)_1)*p_((S_i)_2))*qq_(i,0);
	sum apply(OE,k->p_k*M1_(position(OE,l->l==k),position(OE,l->l==(S_i)_0))
	                           *M2_(position(OE,l->l==k),position(OE,l->l==(S_i)_1))
		                   *M3_(position(OE,l->l==k),position(OE,l->l==(S_i)_2)))
) 
qq=matrix qq;
nonrepeatedTN=time unique(flatten entries qq);
length nonrepeatedTN --34 --> 64-34=30 equalities

qq_(position(S,i->i==(A,C,T)),0)==qq_(position(S,i->i==(A,T,C)),0) --true
pp_(position(S,i->i==(C,C,T)),0)==pp_(position(S,i->i==(T,C,C)),0) --false
qq_(position(S,i->i==(A,C,T)),0)==qq_(position(S,i->i==(G,C,T)),0) --true
qq_(position(S,i->i==(G,C,T)),0)==qq_(position(S,i->i==(G,T,C)),0) --true

qq_(position(S,i->i==(A,C,T)),0)==qq_(position(S,i->i==(G,C,T)),0) --true

qq_(position(S,i->i==(C,G,T)),0)==qq_(position(S,i->i==(T,G,C)),0)--true
qq_(position(S,i->i==(C,A,T)),0)==qq_(position(S,i->i==(T,A,C)),0)--true
qq_(position(S,i->i==(C,G,A)),0)==qq_(position(S,i->i==(C,A,G)),0)--true
qq_(position(S,i->i==(C,C,T)),0)==qq_(position(S,i->i==(T,C,C)),0)--false
qq_(position(S,i->i==(C,T,T)),0)==qq_(position(S,i->i==(T,T,C)),0)--false
qq_(position(S,i->i==(C,C,T)),0)==qq_(position(S,i->i==(T,T,C)),0)--false
p_T*qq_(position(S,i->i==(C,C,T)),0)==p_C*qq_(position(S,i->i==(T,T,C)),0)--false

--b=c=e
EI=sub(qq,flatten toList apply(1..3,i->{c_i=>b_i,e_i=>b_i}));
nonrepeatedEI=time unique(flatten entries EI);
length nonrepeatedEI --17 --> 64-17=47 equalities

--c=e
HKY=sub(qq,flatten toList apply(1..3,i->{e_i=>c_i}));
nonrepeatedHKY=time unique(flatten entries HKY);
length nonrepeatedHKY --34

HKY_(position(S,i->i==(A,C,T)),0)==HKY_(position(S,i->i==(A,T,C)),0) 
HKY_(position(S,i->i==(A,C,T)),0)==HKY_(position(S,i->i==(G,C,T)),0) 


