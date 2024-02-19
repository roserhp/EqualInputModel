-- 5 leaves
-- tree topology: 12|5|34
-- change of basis: H=transpose matrix{{1,1,1,1},{0,1/p_C,0,-1/p_T},{1/(p_A+p_G),-1/(p_C+p_T),1/(p_A+p_G),-1/(p_C+p_T)},{1/p_A,0,-1/p_G,0}}
-- parametrization on the eigenvalues (basis H): 
-- qbar: with identity at the leaves, 
-- pbar: general
restart
K=frac(QQ[p_A,p_C,p_G,p_T]);
Rl=K[l_(6,A),l_(6,C),l_(6,G),l_(6,T),l_(7,A),l_(7,C),l_(7,G),l_(7,T)]
qbar=value get "5leaves_tensor_id_basisH.txt";

OE={A,C,G,T}
S=sort elements (set OE)^**5/splice/splice/splice;

nonMonomial=select(S,i->(length terms qbar_(position(S,j->j==i),0)>1));
length nonMonomial --88
netList nonMonomial
--Example:
qbar_(position(S,j->j==(T,T,T,T,T)),0)

nonZeroEntries=S_(positions(flatten entries qbar,i->i!=0));
monomialNonZeroEntries=select(nonZeroEntries,i->(length terms qbar_(position(S,j->j==i),0)==1));
length monomialNonZeroEntries --225
--Example:
qbar_(position(S,j->j==(A,A,A,T,T)),0)

Rgeneral=K[l_(1,A),l_(1,C),l_(1,G),l_(1,T),l_(2,A),l_(2,C),l_(2,G),l_(2,T),l_(3,A),l_(3,C),l_(3,G),l_(3,T),l_(4,A),l_(4,C),l_(4,G),l_(4,T),l_(5,A),l_(5,C),l_(5,G),l_(5,T),l_(6,A),l_(6,C),l_(6,G),l_(6,T),l_(7,A),l_(7,C),l_(7,G),l_(7,T)]
qbar=sub(qbar,Rgeneral);
pbar=transpose matrix{toList apply(S,i->l_(1,i_0)*l_(2,i_1)*l_(3,i_2)*l_(4,i_3)*l_(5,i_4)*qbar_(position(S,j->j==i),0))};

--Flattening 12|345

s={(A,A),(A,T),(T,A),(T,G),(G,T),(T,T),(G,G),(A,G),(G,A),(C,C),(A,C),(C,A),(C,G),(G,C),(T,C),(C,T)}
sA=toList apply(s,i->append(i,A));
sC=toList apply(s,i->append(i,C));
sG=toList apply(s,i->append(i,G));
sT=toList apply(s,i->append(i,T));
srow=s;
scol=sA|sG|sC|sT;
flattq=mutableMatrix(Rl,4^2,4^3);
for i to 15 do (for j to 63 do
    flattq_(i,j)=qbar_(position(S,k->k==((srow_i)_0,(srow_i)_1,(scol_j)_0,(scol_j)_1,(scol_j)_2)),0);
    );
flattq=matrix flattq;
--sanity check
flattq_(2,3)==qbar_(position(S,j->j==(T,A,T,G,A)),0)
flattq_(2,19)==qbar_(position(S,j->j==(T,A,T,G,G)),0)

flattqPI=sub(flattq,{p_A=>1/2,p_C=>1/3,p_G=>1/8,p_T=>1/24});
rank flattqPI --4
FqA=flattqPI_(toList (0..15))
rank FqA --4
FqG=flattqPI_(toList (16..31))
rank FqG --4
FqC=flattqPI_(toList (32..47))
rank FqC --2
FqT=flattqPI_(toList (48..63))
rank FqT --2
flattqPID=sub(flattqPI,apply(gens Rl,i->i=>1));

--General case
Rgeneral=K[l_(1,A),l_(1,C),l_(1,G),l_(1,T),l_(2,A),l_(2,C),l_(2,G),l_(2,T),l_(3,A),l_(3,C),l_(3,G),l_(3,T),l_(4,A),l_(4,C),l_(4,G),l_(4,T),l_(5,A),l_(5,C),l_(5,G),l_(5,T),l_(6,A),l_(6,C),l_(6,G),l_(6,T),l_(7,A),l_(7,C),l_(7,G),l_(7,T)]
qbar=sub(qbar,Rgeneral);
pbar=transpose matrix{toList apply(S,i->l_(1,i_0)*l_(2,i_1)*l_(3,i_2)*l_(4,i_3)*l_(5,i_4)*qbar_(position(S,j->j==i),0))};

flattp=mutableMatrix(Rgeneral,4^2,4^3);
for i to 15 do (for j to 63 do
    flattp_(i,j)=pbar_(position(S,k->k==((srow_i)_0,(srow_i)_1,(scol_j)_0,(scol_j)_1,(scol_j)_2)),0);
    );
flattp=matrix flattp;
--sanity check
flattp_(2,3)==pbar_(position(S,j->j==(T,A,T,G,A)),0)
flattp_(2,19)==pbar_(position(S,j->j==(T,A,T,G,G)),0)

flattpPI=sub(flattp,{p_A=>1/2,p_C=>1/3,p_G=>1/8,p_T=>1/24});
rank flattpPI --4
FpA=flattpPI_(toList (0..15))
rank FpA --4
rank FpA_{0,1,7,10} --(A,A,A),(A,T,A),(A,G,A),(A,C,A)
FpG=flattpPI_(toList (16..31))
rank FpG --4
FpC=flattpPI_(toList (32..47))
rank FpC --2
FpT=flattpPI_(toList (48..63))
rank FpT --2
rank FpT_{0,1} --(A,A,T),(A,T,T)

flattpPID=sub(flattpPI,apply(gens Rgeneral,i->i=>1));
flattpPID
rank flattpPID --4
FpAid=flattpPID_(toList (0..15))
rank FpAid --4
rank FpAid_{0,1,7,10} --(A,A,A),(A,T,A),(A,G,A),(A,C,A)
FpGid=flattpPID_(toList (16..31))
rank FpGid --4
FpCid=flattpPID_(toList (32..47))
rank FpCid --2
FpTid=flattpPID_(toList (48..63))
rank FpTid --2

------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
--COMPUTATION OF qbar (only needed the first time)
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------

restart
K=frac(QQ[p_A,p_C,p_G,p_T]);
Rl=K[l_(6,A),l_(6,C),l_(6,G),l_(6,T),l_(7,A),l_(7,C),l_(7,G),l_(7,T)]
Hl=transpose(matrix{{1,1,1,1},{0,1/p_C,0,-1/p_T},{1/(p_A+p_G),-1/(p_C+p_T),1/(p_A+p_G),-1/(p_C+p_T)},{1/p_A,0,-1/p_G,0}});
--interior edges
M6=Hl*diagonalMatrix(Rl,4,4,{l_(6,A),l_(6,C),l_(6,G),l_(6,T)})*inverse(Hl)
M7=Hl*diagonalMatrix(Rl,4,4,{l_(7,A),l_(7,C),l_(7,G),l_(7,T)})*inverse(Hl)
--sanity check
inverse(Hl)*M6*Hl
inverse(Hl)*M7*Hl

--identity at the leaves
M1=id_(Rl^4)
M2=M1
M3=M1
M4=M1
M5=M1

OE={A,C,G,T}

S=sort elements (set OE)^**5/splice/splice/splice

qq=mutableMatrix(Rl,4^5,1)

--Topology 12|5|34
for i to (4^5-1) do (
	qq_(i,0)=sum flatten flatten apply(OE,r->apply(OE,s->apply(OE,t->p_r*
	    M1_(position(OE,l->l==r),position(OE,l->l==(S_i)_0))*
	    M2_(position(OE,l->l==r),position(OE,l->l==(S_i)_1))*
	    M3_(position(OE,l->l==t),position(OE,l->l==(S_i)_2))*
	    M4_(position(OE,l->l==t),position(OE,l->l==(S_i)_3))*
	    M5_(position(OE,l->l==s),position(OE,l->l==(S_i)_4))*
	    M6_(position(OE,l->l==r),position(OE,l->l==s))*
	    M7_(position(OE,l->l==s),position(OE,l->l==t)))))
) 

qq=matrix qq;

length select(flatten entries qq,i->i!=0) --64 out of 1024
netList select(flatten entries qq,i->i!=0)

--Orthogonal basis w.r.t. pi-scalar product
H=transpose(matrix{{1,1,1,1},{0,1/p_C,0,-1/p_T},{1/(p_A+p_G),-1/(p_C+p_T),1/(p_A+p_G),-1/(p_C+p_T)},{1/p_A,0,-1/p_G,0}});
H5=(transpose H)**(transpose H)**(transpose H)**(transpose H)**(transpose H);
qbar=time H5*qq;
-- used 421.312 seconds
"5leaves_tensor_id_basisH.txt" << toString qbar << endl << close

