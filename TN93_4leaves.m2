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

