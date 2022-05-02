new;
_dxmiss=0/0;

/*Aquaro, M., Bailey, N. and Pesaran, H. P. (2020). A Quasi Maximum Likelihood Estimation of Spatial Models with Heterogeneous Parameters, Journal of Applied Econometrics, forthcoming */

/*********** Input Data *******************/
/* Dataset N=377, T=159 */
load x_raw[377,159]="hp_raw_new.txt"; /* House Prices: raw data */
load x_res[377,159]="hp_defact_new.txt"; /* House Prices: residuals */

/********** Data Preparation ************/
a_size=0.05;
delta=1/2;

/* Choice of dataset */
x_raw=x_raw';
x_res=x_res';

/* Standardise data */
{x_rawst}=standx(x_raw);
{x_resst}=standx(x_res);

/*********** CD test *************/
{av_c_hp,cd_hp}=CD(x_raw);
{av_c_res,cd_res}=CD(x_res);

/********** Exponent of cross-sectional dependence (2016 version) ************/
{a_hp,se_hp}=atildeall(x_rawst,a_size);

/* ***************************** Output ******************************* */
"";"";
"CD statistic, Exponent of CSD (alpha) and s.e. ";
"House prices";cd_hp~a_hp~se_hp;
"Residual house prices";cd_res;
"";"";

/* ******************************************************** Procedures ******************************************************* */
/* Standardise data */
proc standx(x);
local n,t,m_x,std_x,x_stand;
n=cols(x);
t=rows(x);
m_x=meanc(x)';
std_x=stdc(x)';
x_stand=zeros(t,n);
for i(1,n,1);
x_stand[.,i]=(x[.,i]-m_x[1,i])/std_x[1,i];
endfor;
retp(x_stand);
endp;

/* Pesaran, M.H. (2015). Testing Weak Cross-sectional Dependence in Large Panels. Econometrics Reviews, 34(6-10), 1089-1117 */
/* CD test */
proc (2)=CD(x);
local n,t,corr_m,low_diag,seq_units,pair_corr,avg_p_corr,wdc,reject_wdc;
n=cols(x);
t=rows(x);
corr_m=corrx(x);
low_diag=vech(corr_m);
seq_units=seqa(1,1,rows(low_diag));
pair_corr=rev(sortc(low_diag~seq_units,1));
pair_corr=pair_corr[n+1:rows(pair_corr),1];
avg_p_corr=(2/(n*(n-1)))*sumc(pair_corr);
wdc=sqrt((t)/2)*sqrt(n*(n-1))*avg_p_corr;
if abs(wdc)>1.96;
reject_wdc=1;
else;
reject_wdc=0;
endif;
retp(avg_p_corr,wdc);
endp;

/* Bailey, N., Kapetanios, G. and Pesaran, H. P. (2016) Exponent of cross-sectional Dependence: Estimation and Inference. Journal of Applied Econometrics 31(6), 929-960.

Temporal structure on factors; cross-sectional dependence on errors
{a_tilde,a_thrtilde,omega_tild,omega_thrtild}=atildeall(x,a_size)

Input: x (TxN), a_size
N: cross section dimension
T: time series dimension
a_size: scalar for significance level of thresholding: set to eg. 0.01, 0.05, 0.10 

Output:
Bias corrected estimate of alpha and respective s.e.
 */
/* a_tilde estimates: temporal structure on factors; cross-sectional dependence on errors */
proc (2)=atildeall(x,a_size);
local n,t,p,z,ln_z,x_bar1,x_bar1_c,std_x_bar1,pc,ln_var,x_bar,m_x_bar,x_bar_stand,i,x_bar_2m,m_x_bar_2m,x_bar_2m_st,x_bar_2m_st_lag,
v_all,v_1,dv,rhs,b,s_b,e_nw,sse_nw,sig2_nw,v_f_2,c_avg,m_c_avg,c_avg_2,index,e,e_bar,m_e_bar,e_bar_stand,e_bar_2m,v_e_bar,s_hat,a_dot,
a_tilde,rhsx,coefx,residx,ssex,sdx,dfx,t_test,p_n,theta,size,x_str,x_str1,xstr_bar,musqr_thr,a_thrtilde,ggg_o,c_avg_sel,m_c_avg_sel,
frasel,s_frasel,ggg_t,c_avg_selt,m_c_avg_selt,fraselt,s_fraselt,j,s_ttest,s_size,omega_tild;

n=cols(x);
t=rows(x);
p=ceil(t^(1/3));
z=seqa(1,1,n); /* cross sectional trend */
ln_z=ln(z);
x_bar1=meanc(x');
x_bar1_c=x_bar1;
std_x_bar1=stdc(x_bar1_c);
ln_var=ln(std_x_bar1^2);
x_bar=x_bar1_c./std_x_bar1; /* standardise the cross-sectional avgs */
m_x_bar=meanc(x_bar);
x_bar_stand=zeros(t,1);
i=1;
do while i<=t;
x_bar_stand[i,.]=x_bar[i,.]-m_x_bar;
i=i+1;
endo;
/*Newey-West method*/
x_bar_2m=x_bar_stand^2;
m_x_bar_2m=meanc(x_bar_2m);
x_bar_2m_st=zeros(t,1);
for i(1,t,1);
x_bar_2m_st[i,.]=x_bar_2m[i,.]-m_x_bar_2m;
endfor;
x_bar_2m_st_lag=zeros(rows(x_bar_2m_st),p);
for i(1,p,1);
x_bar_2m_st_lag[.,i]=lagn(x_bar_2m_st,i);
endfor;
v_all=x_bar_2m_st~x_bar_2m_st_lag;
v_1=v_all[p+1:t,.];
dv=v_1[.,1];
rhs=v_1[.,2:p+1];
b=inv(rhs'*rhs)*rhs'*dv;
s_b=sumc(b);
e_nw=dv-rhs*b;
sse_nw=e_nw'*e_nw;
sig2_nw=sse_nw/(t-cols(rhs));
v_f_2=sig2_nw/(1-s_b)^2;

c_avg=inv(x_bar'*x_bar)*x_bar'*x; /* OLS estimate standardised cross-sectional coefficients*/
m_c_avg=meanc(c_avg');
{pc}=getpc(x,4);
c_avg_2=inv(pc'*pc)*pc'*x; /* OLS estimate non standardised cross-sectional coefficients*/
index=rev(sortc(abs(c_avg_2')~z,1));
e=x-pc*c_avg_2; /* calculate residuals from non standardised cross-sectionals regression */
e_bar=meanc(e');
m_e_bar=meanc(e_bar);
e_bar_stand=zeros(t,1);
i=1;
do while i<=t;
e_bar_stand[i,.]=e_bar[i,.]-m_e_bar;
i=i+1;
endo;
e_bar_2m=e_bar_stand^2;
v_e_bar=meanc(e_bar_2m);
s_hat=n*v_e_bar;

a_dot=1+(1/2)*(ln_var/ln(n));
a_tilde=a_dot-(1/2)*(s_hat/(n*ln(n)*std_x_bar1^2));

rhsx=ones(t,1)~x_bar1;
coefx=inv(rhsx'*rhsx)*rhsx'*x;
residx=x-rhsx*coefx;
ssex=residx'*residx/(t-cols(rhsx));
sdx=sqrt(diag(ssex*inv(x_bar1'*x_bar1)));
dfx=t-cols(rhsx);
t_test=coefx[2,.]'./sdx[.,1];
size=zeros(n,1);
x_str=zeros(t,n);
s_ttest=rev(sortc(abs(t_test)~z,1));
j=1;
do while  j<=cols(coefx);
p_n=a_size/(n-j+1);
theta=cdfni(1-p_n/2);
  if abs(s_ttest[j,1])>=theta;
   size[j,1]=1;
  else;
   size[j,1]=0;
  endif;
j=j+1;
endo;
s_size=sortc(size~s_ttest[.,2],2);
x_str=s_size[.,1]'.*x;
x_str=x_str';
x_str1=delif(x_str,x_str[.,1].==0);
if x_str1==miss(1,1);
   musqr_thr=1;
  else;
   x_str1=x_str1';
   xstr_bar=meanc(x_str1');
   musqr_thr=meanc((xstr_bar-meanc(xstr_bar))^2);
endif;
a_thrtilde=a_tilde-(1/2)*(ln(musqr_thr)/ln(n));

ggg_o=round(n^(a_tilde));
    if ggg_o>=n;
       ggg_o=n;
    elseif ggg_o<1;
        ggg_o=1;
    else;
       ggg_o=ggg_o;
    endif;
c_avg_sel=rev(sortc(c_avg'~abs(c_avg'),2));
c_avg_sel=c_avg_sel[1:ggg_o,1];
m_c_avg_sel=meanc(c_avg_sel);
frasel=zeros(1,ggg_o);
for i(1,ggg_o,1);
  frasel[.,i]=(c_avg_sel[i,1]-m_c_avg_sel)^2;
endfor;
s_frasel=sumc(frasel');

ggg_t=round(n^(a_thrtilde));
    if ggg_t>=n;
       ggg_t=n;
    elseif ggg_t<1;
        ggg_t=1;
    else;
       ggg_t=ggg_t;
    endif;
c_avg_selt=rev(sortc(c_avg'~abs(c_avg'),2));
c_avg_selt=c_avg_selt[1:ggg_t,1];
m_c_avg_selt=meanc(c_avg_selt);
fraselt=zeros(1,ggg_t);
for i(1,ggg_t,1);
  fraselt[.,i]=(c_avg_selt[i,1]-m_c_avg_selt)^2;
endfor;
s_fraselt=sumc(fraselt');
omega_tild=((1/t)*(v_f_2)+(4/n)*(n^(1-a_thrtilde)*s_fraselt/(ggg_t-1)))^(1/2)/(2*ln(n)); /* 90% CI: 1.65, 95% CI: 1.96 */ 
retp(a_thrtilde,omega_tild); /*output*/
endp;

proc getpc(ax,j1);
local j,x,n,t,first_j,eigvals,eigvecs,evals,evecs,pc;
x=ax;
n=cols(ax);
t=rows(ax);
j=j1; /* or 'factor' if more than the maximum pc is needed */
if n<t;
first_j=seqa(n,-1,j);
{eigvals,eigvecs}=eigrs2(x'x); /*sorted in ascending order*/
evals=submat(eigvals,first_j,1); /* pick j largest eigenvalues */
evecs=submat(eigvecs,0,first_j);
pc=x*evecs;
elseif n>=t;
first_j=seqa(t,-1,j);
{eigvals,eigvecs}=eigrs2(x*x'); /*sorted in ascending order*/
evals=submat(eigvals,first_j,1); /* pick j largest eigenvalues */
evecs=submat(eigvecs,0,first_j);
pc=evecs; /* eigenvectors of xx' = PC of x */
endif;
retp(pc);
endp;
