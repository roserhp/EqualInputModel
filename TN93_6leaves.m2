restart
K=frac(QQ[p_A,p_C,p_G,p_T]);
R=K[b_1,c_1,e_1,b_2,c_2,e_2,b_3,c_3,e_3,b_4,c_4,e_4,b_5,c_5,e_5,b_6,c_6,e_6,b_7,c_7,e_7,b_8,c_8,e_8,b_9,c_9,e_9];
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
i=8
M8=sub(M1,{b_1=>b_i,c_1=>c_i,e_1=>e_i})
i=9
M9=sub(M1,{b_1=>b_i,c_1=>c_i,e_1=>e_i})

--identidad en las hojas
M1=id_(R^4)
M2=M1
M3=M1
M4=M1
M5=M1
M6=M1

{M1,M2,M3,M4,M5,M6,M7,M8,M9}

S=sort elements (set OE)^**6/splice/splice/splice/splice

pp=mutableMatrix(R,4^6,1)

--Flattening 12|34|56
time for i to (4^6-1) do (
	pp_(i,0)=sum flatten flatten flatten apply(OE,r->apply(OE,s->apply(OE,t->apply(OE,u->p_r*
	    M1_(position(OE,l->l==r),position(OE,l->l==(S_i)_0))*
	    M2_(position(OE,l->l==r),position(OE,l->l==(S_i)_1))*
	    M3_(position(OE,l->l==t),position(OE,l->l==(S_i)_2))*
	    M4_(position(OE,l->l==t),position(OE,l->l==(S_i)_3))*
	    M5_(position(OE,l->l==u),position(OE,l->l==(S_i)_4))*
	    M6_(position(OE,l->l==u),position(OE,l->l==(S_i)_5))*
	    M7_(position(OE,l->l==r),position(OE,l->l==s))*
	    M8_(position(OE,l->l==s),position(OE,l->l==t))*
	    M9_(position(OE,l->l==s),position(OE,l->l==u))))))
) 
-- used 187.578 seconds
pp=matrix pp;


-- NEW BASIS
H=transpose(matrix{{4,4,4,4},{0,1/p_C,0,-1/p_T},{1/(p_A+p_G),-1/(p_C+p_T),1/(p_A+p_G),-1/(p_C+p_T)},{1/p_A,0,-1/p_G,0}});

H6=time (transpose H)**(transpose H)**(transpose H)**(transpose H)**(transpose H)**(transpose H);
 -- used 97.9687 seconds

pbar=time H6*pp;
-- used 554.766 seconds

length select(flatten entries pbar,i->i==0) --2918
zeroEntries=S_(positions(flatten entries pbar,i->i==0));


"6leavesBasisHzeroEntries3Cherries.txt" << netList zeroEntries << endl << close

nonZeroEntries=S_(positions(flatten entries pbar,i->i!=0));
length unique nonZeroEntries --1178
length nonZeroEntries --1178

"6leavesBasisHnonZeroEntries3Cherries.txt" << netList nonZeroEntries << endl << close

--Checking conjecture of non-zeroes 
T2=nonZeroEntries_(positions(nonZeroEntries,i->number(i,j->j==T)>=2 and member(C,i)==false))
C2=nonZeroEntries_(positions(nonZeroEntries,i->number(i,j->j==C)>=2 and member(T,i)==false))

T2bar=S_(positions(S,i->number(i,j->j==T)>=2 and member(C,i)==false))
C2bar=S_(positions(S,i->number(i,j->j==C)>=2 and member(T,i)==false))

C2==C2bar --true: CC..... with no additional T's
T2==T2bar --true: TT..... with no additional C's

length C2 --473
length T2 --473

G2bar=S_(positions(S,i->member(C,i)==false and member(T,i)==false and number(i,j->j==G)>=2))
G2=nonZeroEntries_(positions(nonZeroEntries,i->member(C,i)==false and member(T,i)==false and number(i,j->j==G)>=2))

G2bar==G2 --true: GG..... with additional G's or A's

length G2 --57

member((A,A,A,A,A,A),nonZeroEntries) --true: AAAAAA

1178-473-473-57-1 --174

CCTTbar=S_(positions(S,i->number(i,j->j==C)>=2 and number(i,j->j==T)>=2))
length CCTTbar --650
CCTT=nonZeroEntries_(positions(nonZeroEntries,i->number(i,j->j==C)>=2 and number(i,j->j==T)>=2))
length CCTT --174

a0=select(CCTT,i->(i_0==C and i_1==C) or (i_0==T and i_1==T) or (i_2==C and i_3==C) or (i_2==T and i_3==T)) 
a1=select(CCTT,i->(i_0==C and i_1==C) or (i_0==T and i_1==T) or (i_2==C and i_3==C) or (i_2==T and i_3==T) or (i_4==C and i_5==C) or (i_4==T and i_5==T)) 
a2=select(CCTTbar,i->(i_0==C and i_1==C) or (i_0==T and i_1==T) or (i_2==C and i_3==C) or (i_2==T and i_3==T) or (i_4==C and i_5==C) or (i_4==T and i_5==T)) 

a0==a1 --false
length a0 --142
length a1 --174
length a2 --306

netList sort elements (set a2-set a1) --cases where the some cherry is CT
netList sort elements set a1


-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------
pp2=mutableMatrix(R,4^6,1)

--Flattening 12|34|5|6
time for i to (4^6-1) do (
	pp2_(i,0)=sum flatten flatten flatten apply(OE,r->apply(OE,s->apply(OE,t->apply(OE,u->p_r*
	    M1_(position(OE,l->l==r),position(OE,l->l==(S_i)_0))*
	    M2_(position(OE,l->l==r),position(OE,l->l==(S_i)_1))*
	    M3_(position(OE,l->l==t),position(OE,l->l==(S_i)_2))*
	    M4_(position(OE,l->l==t),position(OE,l->l==(S_i)_3))*
	    M5_(position(OE,l->l==s),position(OE,l->l==(S_i)_4))*
	    M6_(position(OE,l->l==u),position(OE,l->l==(S_i)_5))*
	    M7_(position(OE,l->l==r),position(OE,l->l==s))*
	    M8_(position(OE,l->l==s),position(OE,l->l==u))*
	    M9_(position(OE,l->l==u),position(OE,l->l==t))))))
) 
-- used 201.656 seconds
pp2=matrix pp2;

pbar2=time H6*pp2;
-- used 513.047 seconds

length select(flatten entries pbar2,i->i==0) --2916
zeroEntries2=S_(positions(flatten entries pbar2,i->i==0));


"6leavesBasisHzeroEntries2Cherries.txt" << netList zeroEntries2 << endl << close

nonZeroEntries2=S_(positions(flatten entries pbar2,i->i!=0));
length unique nonZeroEntries2 --1180
length nonZeroEntries2 --1180

"6leavesBasisHnonZeroEntries2Cherries.txt" << netList nonZeroEntries2 << endl << close

--Checking conjecture of non-zeroes 
T22=nonZeroEntries2_(positions(nonZeroEntries2,i->number(i,j->j==T)>=2 and member(C,i)==false))
C22=nonZeroEntries2_(positions(nonZeroEntries2,i->number(i,j->j==C)>=2 and member(T,i)==false))

C22==C2bar --true: CC..... with no additional T's
T22==T2bar --true: TT..... with no additional C's

length C22 --473
length T22 --473

G22=nonZeroEntries2_(positions(nonZeroEntries2,i->member(C,i)==false and member(T,i)==false and number(i,j->j==G)>=2))

G2bar==G22 --true: GG..... with additional G's or A's

length G22 --57

member((A,A,A,A,A,A),nonZeroEntries) --true: AAAAAA

length nonZeroEntries2-length C22-length T22-length G22-1 --176

CCTT2=nonZeroEntries2_(positions(nonZeroEntries2,i->number(i,j->j==C)>=2 and number(i,j->j==T)>=2))
length CCTT2 --176

b0=select(CCTT2,i->(i_0==C and i_1==C) or (i_0==T and i_1==T) or (i_2==C and i_3==C) or (i_2==T and i_3==T)) 
b1=select(CCTT2,i->(i_0==C and i_1==C) or (i_0==T and i_1==T) or (i_2==C and i_3==C) or (i_2==T and i_3==T) or (i_4==C and i_5==C) or (i_4==T and i_5==T)) 
b2=select(CCTTbar,i->(i_0==C and i_1==C) or (i_0==T and i_1==T) or (i_2==C and i_3==C) or (i_2==T and i_3==T)) 

b0==b1 --true

length b1 --144 --NOT ENOUGH ELEMENTS!!! 
length b2 --234

c1=sort elements (set nonZeroEntries2-set C22-set T22-set G22-set {(A,A,A,A,A,A)}-set b1)
length c1+length b1==length CCTT2 --true
netList c1
d1=select(CCTT2,i->(i_4==C and i_5==T) or (i_4==T and i_5==C)) 
c1==d1
netList sort elements (set d1-set c1)

b11=select(CCTT2,i->(i_0==C and i_1==C) or (i_0==T and i_1==T) or (i_2==C and i_3==C) or (i_2==T and i_3==T) or (i_4==C and i_5==T) or (i_4==T and i_5==C));
length b11 --176
b22=select(CCTTbar,i->(i_0==C and i_1==C) or (i_0==T and i_1==T) or (i_2==C and i_3==C) or (i_2==T and i_3==T) or (i_4==C and i_5==T) or (i_4==T and i_5==C));
length b22 --402

netList sort elements (set b22-set b11) 

b3=select(CCTTbar,i->(i_4==C and i_5==T) or (i_4==T and i_5==C));
length b3
b4=select(CCTT2,i->(i_4==C and i_5==T) or (i_4==T and i_5==C));
length b4

netList sort elements (set b3-set b4) 

netList sort elements (set b3-set b2) 

length elements(set b3-set b4) --170
length elements(set b3-set b2) --168
length elements((set b3-set b4)-(set b3-set b2)) --34
