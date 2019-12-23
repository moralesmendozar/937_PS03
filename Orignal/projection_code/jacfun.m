function symJac = jacfun(R,q,cF,cG,kGpol,muFplus,muFminus,lamFplus,lamFminus,lamGplus,lamGminus,a,alower,kF,levF,cF_next1,cF_next2,cF_next3,cF_next4,cF_next5,cG_next1,cG_next2,cG_next3,cG_next4,cG_next5,q_next1,q_next2,q_next3,q_next4,q_next5,kFpol_next1,kFpol_next2,kFpol_next3,kFpol_next4,kFpol_next5,prnext1,prnext2,prnext3,prnext4,prnext5)
%JACFUN
%    SYMJAC = JACFUN(R,Q,CF,CG,KGPOL,MUFPLUS,MUFMINUS,LAMFPLUS,LAMFMINUS,LAMGPLUS,LAMGMINUS,A,ALOWER,KF,LEVF,CF_NEXT1,CF_NEXT2,CF_NEXT3,CF_NEXT4,CF_NEXT5,CG_NEXT1,CG_NEXT2,CG_NEXT3,CG_NEXT4,CG_NEXT5,Q_NEXT1,Q_NEXT2,Q_NEXT3,Q_NEXT4,Q_NEXT5,KFPOL_NEXT1,KFPOL_NEXT2,KFPOL_NEXT3,KFPOL_NEXT4,KFPOL_NEXT5,PRNEXT1,PRNEXT2,PRNEXT3,PRNEXT4,PRNEXT5)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    20-Dec-2019 17:48:08

t2 = cG_next1.^2;
t3 = cG_next2.^2;
t4 = cG_next3.^2;
t5 = cG_next4.^2;
t6 = cG_next5.^2;
t7 = 1.0./cF_next1.^2;
t8 = 1.0./cF_next2.^2;
t9 = 1.0./cF_next3.^2;
t10 = 1.0./cF_next4.^2;
t11 = 1.0./cF_next5.^2;
t17 = 1.0./kGpol.^(3.0./1.0e+1);
t12 = 1.0./t2;
t13 = 1.0./t3;
t14 = 1.0./t4;
t15 = 1.0./t5;
t16 = 1.0./t6;
t18 = prnext1.*t7.*(1.9e+1./2.0e+1);
t19 = prnext2.*t8.*(1.9e+1./2.0e+1);
t20 = prnext3.*t9.*(1.9e+1./2.0e+1);
t21 = prnext4.*t10.*(1.9e+1./2.0e+1);
t22 = prnext5.*t11.*(1.9e+1./2.0e+1);
t23 = t18+t19+t20+t21+t22;
symJac = reshape([-1.0./R.^2,0.0,0.0,-cF.^2.*t23,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,muFplus.*(-9.0./1.0e+1)+1.0,kF+kGpol./1.0e+1-kF.*levF-1.0./1.0e+1,0.0,0.0,0.0,0.0,0.0,R.*cF.*t23.*-2.0,cF.*(t18.*(q_next1+9.801333289067057e-1)+t19.*(q_next2+9.8998383262698e-1)+t21.*(q_next4+1.009982832650868)+t22.*(q_next5+1.020133328871151)+t20.*(q_next3+9.999333355555061e-1)).*-2.0,-1.0,0.0,-1.0,cG.*(prnext1.*t12.*(4.9e+1./5.0e+1)+prnext2.*t13.*(4.9e+1./5.0e+1)+prnext3.*t14.*(4.9e+1./5.0e+1)+prnext4.*t15.*(4.9e+1./5.0e+1)+prnext5.*t16.*(4.9e+1./5.0e+1)).*-2.0,cG.*(prnext4.*t15.*(q_next4+t17.*6.36289184570047e-1).*(4.9e+1./5.0e+1)+prnext3.*t14.*(q_next3+t17.*6.299580013999688e-1).*(4.9e+1./5.0e+1)+prnext5.*t16.*(q_next5+t17.*6.426839971888252e-1).*(4.9e+1./5.0e+1)+prnext2.*t13.*(q_next2+t17.*6.236898145549974e-1).*(4.9e+1./5.0e+1)+prnext1.*t12.*(q_next1+t17.*6.174839972112245e-1).*(4.9e+1./5.0e+1)).*-2.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,cG.^2.*1.0./kGpol.^(1.3e+1./1.0e+1).*t12.*t13.*t14.*t15.*t16.*(prnext1.*t3.*t4.*t5.*t6.*5.56180139949543e+15+prnext2.*t2.*t4.*t5.*t6.*5.61769843284932e+15+prnext3.*t2.*t3.*t5.*t6.*5.674157240727924e+15+prnext4.*t2.*t3.*t4.*t6.*5.731183469059081e+15+prnext5.*t2.*t3.*t4.*t5.*5.788782820513148e+15).*3.26405569239796e-17,1.0,0.0,0.0,q./1.0e+1,-1.0,0.0,0.0,0.0,0.0,-1.0,q.*(-9.0./1.0e+1),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0],[8,11]);
