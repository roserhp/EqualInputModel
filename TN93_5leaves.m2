restart
K=frac(QQ[p_A,p_C,p_G,p_T]);
R=K[b_1,c_1,e_1,b_2,c_2,e_2,b_3,c_3,e_3,b_4,c_4,e_4,b_5,c_5,e_5,b_6,c_6,e_6,b_7,c_7,e_7];
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
i=6
M6=sub(M1,{b_1=>b_i,c_1=>c_i,e_1=>e_i})
i=7
M7=sub(M1,{b_1=>b_i,c_1=>c_i,e_1=>e_i})

--identidad en las hojas
M1=id_(R^4)
M2=M1
M3=M1
M4=M1
M5=M1

{M1,M2,M3,M4,M5,M6,M7}

S=sort elements (set OE)^**5/splice/splice/splice

pp=mutableMatrix(R,4^5,1)

--Flattening 12|5|34
for i to (4^5-1) do (
	pp_(i,0)=sum flatten flatten apply(OE,r->apply(OE,s->apply(OE,t->p_r*
	    M1_(position(OE,l->l==r),position(OE,l->l==(S_i)_0))*
	    M2_(position(OE,l->l==r),position(OE,l->l==(S_i)_1))*
	    M3_(position(OE,l->l==t),position(OE,l->l==(S_i)_2))*
	    M4_(position(OE,l->l==t),position(OE,l->l==(S_i)_3))*
	    M5_(position(OE,l->l==s),position(OE,l->l==(S_i)_4))*
	    M6_(position(OE,l->l==r),position(OE,l->l==s))*
	    M7_(position(OE,l->l==s),position(OE,l->l==t)))))
) 

pp=matrix pp;

-- NEW BASIS
H=transpose(matrix{{4,4,4,4},{0,1/p_C,0,-1/p_T},{1/(p_A+p_G),-1/(p_C+p_T),1/(p_A+p_G),-1/(p_C+p_T)},{1/p_A,0,-1/p_G,0}});
H5=(transpose H)**(transpose H)**(transpose H)**(transpose H)**(transpose H);
pbar=time H5*pp;
--33 sec with identity at the leaves

length select(flatten entries pbar,i->i==0) 
zeroEntries=S_(positions(flatten entries pbar,i->i==0))

"5leavesBasisHzeroEntries.txt" << netList zeroEntries << endl << close

nonZeroEntries=S_(positions(flatten entries pbar,i->i!=0))
length unique nonZeroEntries --313
length nonZeroEntries --313

"5leavesBasisHnonZeroEntries.txt" << netList nonZeroEntries << endl << close

--Checking conjecture of non-zeroes 
T2=nonZeroEntries_(positions(nonZeroEntries,i->number(i,j->j==T)>=2 and member(C,i)==false))
C2=nonZeroEntries_(positions(nonZeroEntries,i->number(i,j->j==C)>=2 and member(T,i)==false))

T2bar=S_(positions(S,i->number(i,j->j==T)>=2 and member(C,i)==false))
C2bar=S_(positions(S,i->number(i,j->j==C)>=2 and member(T,i)==false))

C2==C2bar --true: CC..... with no additional T's
T2==T2bar --true: TT..... with no additional C's

length C2 --131
length T2 --131

G2bar=S_(positions(S,i->member(C,i)==false and member(T,i)==false and number(i,j->j==G)>=2))
G2=nonZeroEntries_(positions(nonZeroEntries,i->member(C,i)==false and member(T,i)==false and number(i,j->j==G)>=2))

G2bar==G2 --true: GG..... with additional G's or A's

length G2 --26

member((A,A,A,A,A),nonZeroEntries) --true: AAAAA

313-131-131-26-1 --24

CCTTbar=S_(positions(S,i->number(i,j->j==C)>=2 and number(i,j->j==T)>=2))
length CCTTbar --80
netList CCTTbar
CCTT=nonZeroEntries_(positions(nonZeroEntries,i->number(i,j->j==C)>=2 and number(i,j->j==T)>=2))
length CCTT --24
netList CCTT
-- CCTT_,TTCC_ (2*4=8) is non-zero for this model
-- Check that all those entries are in the non-zero list:
a0=select(CCTT,i->(i_0==C and i_1==C and i_2==T and i_3==T) or (i_0==T and i_1==T and i_2==C and i_3==C)) --true
a00=select(CCTTbar,i->(i_0==C and i_1==C and i_2==T and i_3==T) or (i_0==T and i_1==T and i_2==C and i_3==C)) --true

a0==a00 --true

a1=select(CCTT,i->(i_0==C and i_1==C) or (i_0==T and i_1==T) or (i_2==C and i_3==C) or (i_2==T and i_3==T)) --true
a2=select(CCTTbar,i->(i_0==C and i_1==C) or (i_0==T and i_1==T) or (i_2==C and i_3==C) or (i_2==T and i_3==T)) --true

length a1 
length a2

netList sort elements (set a2-set a1) --cases where the second cherry is CT
netList sort elements set a1


--Topologic invariants for 12|5|34
-- CTCT_,CTTC_,TCCT_,TCTC_ (4*4=16) is zero for this model
CCTT0=zeroEntries_(positions(zeroEntries,i->number(i,j->j==C)>=2 and number(i,j->j==T)>=2))
length CCTT0 --56
netList CCTT0
-- for all those elements, there exists a topology for which they are non-zero
sort elements set CCTTbar == sort elements (set CCTT0 + set CCTT)


-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------
--ORIGINAL BASIS
--Repeated entries?
nonrepeatedpp=time unique(flatten entries pp);
length nonrepeatedpp --53 --> 1024-53 equalities

pp_(position(S,i->i==(A,C,G,T)),0)==pp_(position(S,i->i==(A,T,G,C)),0)--true
pp_(position(S,i->i==(C,G,T,A)),0)==pp_(position(S,i->i==(T,G,C,A)),0)--true
pp_(position(S,i->i==(C,A,G,T)),0)==pp_(position(S,i->i==(T,A,G,C)),0)--true

pp_(position(S,i->i==(C,C,G,T)),0)==pp_(position(S,i->i==(T,C,G,C)),0)--false
pp_(position(S,i->i==(C,G,G,T)),0)==pp_(position(S,i->i==(T,G,G,C)),0)--true
pp_(position(S,i->i==(C,T,G,T)),0)==pp_(position(S,i->i==(T,T,G,C)),0)--false

p_G*pp_(position(S,i->i==(A,A,C,T)),0)==p_A*pp_(position(S,i->i==(A,G,C,T)),0)--false

--not updated to 5 leaves
equalEntries=select(unique for i to 255 list S_(positions(flatten entries pp,l->l==pp_(i,0))),k->length k>1);
length equalEntries
netList equalEntries
"4leavesOriginalBasisEqualEntries.txt" << netList equalEntries << endl << close
