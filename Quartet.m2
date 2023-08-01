restart
--Ring definition for QUARTETS with identity at the leaves and inner egde matrix
-- on p1..p4 and EIGENVALUES l1..l4
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(5,1),l_(5,2),l_(5,3),l_(5,4)]
gens R

--Diagonalizing change of basis
H=transpose(matrix{{1,1,1,1},{1/(p_1+p_2),1/(p_1+p_2),-1/(p_3+p_4),-1/(p_3+p_4)},
	{0,0,1/p_3,-1/p_4},{1/p_1,-1/p_2,0,0}})

--sanity check
inverse (transpose H) --A in paper up to p1+p2+p3+p4=1
M=H*diagonalMatrix(R,4,4,toList(l_(5,1)..l_(5,4)))*inverse(H)

--Identity at the leaves
M1=id_(R^4)
M2=M1
M3=M1
M4=M1

st=toList(1..4)
S=sort elements (set st)^**4/splice/splice
netList S

--Tensor for configuration 12|34
qq=mutableMatrix(R,256,1)
for i to 255 do (
	qq_(i,0)=sum flatten toList apply(st,k->apply(st,kk->p_k*M1_(position(st,l->l==k),position(st,l->l==(S_i)_0))*M2_(position(st,l->l==k),position(st,l->l==(S_i)_1))*M3_(position(st,l->l==kk),position(st,l->l==(S_i)_2))*M4_(position(st,l->l==kk),position(st,l->l==(S_i)_3))*M_(position(st,l->l==k),position(st,l->l==kk))))
) 
qq=matrix qq;
netList (flatten entries qq)

--Tensor with change of basis
H4=(transpose H)**(transpose H)**(transpose H)**(transpose H);
qbar=time H4*qq;
-- used 0.976217 seconds
"4leaves_tensor_id_FinalBasis.txt" << toString qbar << endl << close

--Different types of entries in the tensor in the new basis
nonMonomial=select(S,i->(length terms qbar_(position(S,j->j==i),0)>1))
netList nonMonomial
length nonMonomial
nonZeroEntries=S_(positions(flatten entries qbar,i->i!=0))
length nonZeroEntries 
monomialNonZeroEntries=select(nonZeroEntries,i->(length terms qbar_(position(S,j->j==i),0)==1))
length monomialNonZeroEntries
netList monomialNonZeroEntries

--Flattening 12|34 for identity at the leaves
flattq=mutableMatrix(R,16,16)
for i to 15 do (for j to 15 do 
    (if(j%4==j//4 and i%4==i//4) then flattq_(i,j)=p_(st_(i%4))*M_(i%4,j%4));      
    );
flattq=matrix(flattq);
flattq
--Change of basis of the flattening
flattQ=time (transpose(H)**transpose(H))*flattq*(H**H); 
-- used 0.0976014 seconds
--Quasi-block form (reordering according to s in next paragraph of code)
blockQ=flattQ_{0,3,12,7,13,15,5,1,4,10,2,8,6,9,11,14}^{0,3,12,7,13,15,5,1,4,10,2,8,6,9,11,14}
rank blockQ_{1,2,3,4}--1
rank blockQ_{14,15}--0
rank blockQ_{10,11,12,13}--1
rank blockQ_{7,8}--1

--Ring definition for flattening 12|34 general
Rgeneral=K[l_(1,1)..l_(5,4)]
s={(1,1),(1,4),(4,1),(2,4),(4,2),(4,4),(2,2),(1,2),(2,1),(3,3),(1,3),(3,1),(2,3),(3,2),(3,4),(4,3)}
-- Moving from identity at the leaves to general case
flattQ=sub(blockQ,Rgeneral);
flattP=matrix toList apply(0..15,i->toList apply(0..15,j->l_(1,(s_i)_0)*l_(2,(s_i)_1)*l_(3,(s_j)_0)*l_(4,(s_j)_1)*flattQ_(i,j)));
--rk 0
rank flattP_{position(s,i->i==(3,4)),position(s,i->i==(4,3))}--0
--rk 1
rank flattP_{position(s,i->i==(1,4)),position(s,i->i==(4,1)),position(s,i->i==(2,4)),position(s,i->i==(4,3))}--1
rank flattP_{position(s,i->i==(1,3)),position(s,i->i==(3,1)),position(s,i->i==(3,2)),position(s,i->i==(2,3))}--1
rank flattP_{position(s,i->i==(1,2)),position(s,i->i==(2,1))}--1
--rk 2
rank flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(2,2))}--2
--rk 3
rank flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,4)),position(s,i->i==(4,4))}--3
rank flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,3)),position(s,i->i==(3,3))}--3
--rk 4: all
--column generators:
rank flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,3)),position(s,i->i==(1,4))}--4
--note that it forms a diagonal submatrix with non-zero entries in the diagonal as long as pi is generic and all eigenvalues are non-zero
flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,3)),position(s,i->i==(1,4))}^{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,3)),position(s,i->i==(1,4))}

--Tensor in general case
qbar=sub(qbar,Rgeneral);
pbar=transpose matrix{toList apply(S,i->l_(1,i_0)*l_(2,i_1)*l_(3,i_2)*l_(4,i_3)*qbar_(position(S,j->j==i),0))};
"4leaves_tensor_FinalBasis.txt" << toString pbar << endl << close

--------------------------------------------------------------------
--------------------------------------------------------------------
-------------------------------------------------------------------
-- Alternative computation of flattenings
--------------------------------------------------------------------
--------------------------------------------------------------------
-------------------------------------------------------------------

restart

K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(5,1),l_(5,2),l_(5,3),l_(5,4)]
qbar=value get "4leaves_tensor_id_FinalBasis.txt";

S=sort elements (set {1,2,3,4})^**4/splice/splice;
--Flattening 12|34
s={(1,1),(1,4),(4,1),(2,4),(4,2),(4,4),(2,2),(1,2),(2,1),(3,3),(1,3),(3,1),(2,3),(3,2),(3,4),(4,3)}
flattQ=mutableMatrix(R,4^2,4^2);
for i to 15 do (for j to 15 do
    flattQ_(i,j)=qbar_(position(S,k->k==((s_i)_0,(s_i)_1,(s_j)_0,(s_j)_1)),0);
    );
flattQ=matrix flattQ;

Rgeneral=K[l_(1,1)..l_(5,4)]
qbar=sub(qbar,Rgeneral);
pbar=transpose matrix{toList apply(S,i->l_(1,i_0)*l_(2,i_1)*l_(3,i_2)*l_(4,i_3)*qbar_(position(S,j->j==i),0))};
flattQ=sub(flattQ,Rgeneral);
flattP=matrix toList apply(0..15,i->toList apply(0..15,j->l_(1,(s_i)_0)*l_(2,(s_i)_1)*l_(3,(s_j)_0)*l_(4,(s_j)_1)*flattQ_(i,j)));

testpbar=value get "4leaves_tensor_FinalBasis.txt";
testflattP=mutableMatrix(Rgeneral,16,16);
for i to 15 do (for j to 15 do
    testflattP_(i,j)=testpbar_(position(S,k->k==((s_i)_0,(s_i)_1,(s_j)_0,(s_j)_1)),0);
    );
testflattP=matrix testflattP;

testflattP==flattP --true

