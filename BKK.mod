/* This model follows Backus, Kehoe and Kydland (1993). It is a two-country,
two-intermediate good economy with TFP shocks at each country.
A typical IBC setup.
*/


@#define countries = [ "h", "f" ]

@#for co in countries
var C_@{co} N_@{co} K_@{co} X_@{co} A_@{co} B_@{co} Y_@{co} nx_@{co} ir_@{co} rx_@{co} 
qa_@{co} qb_@{co} p_@{co} z_@{co} c_@{co} n_@{co} k_@{co} x_@{co} a_@{co} b_@{co} y_@{co};

varexo eta_@{co};

parameters beta_@{co} mu_@{co} gamma_@{co} theta_@{co} delta_@{co} sigma_@{co} omega_@{co};
@#endfor

parameters rho_h_h rho_h_f rho_f_h rho_f_f;


/* main model starts */
model;
/* intratemporal choice between C_h and N_h */
mu_h*C_h^(mu_h-1)*(C_h^mu_h*(1-N_h)^(1-mu_h))^(gamma_h-1)*omega_h*A_h^(sigma_h-1)*(omega_h*A_h^(sigma_h)+(1-omega_h)*B_h^(sigma_h))^((sigma_h^-1)-1)*exp(z_h)*K_h(-1)^theta_h*(1-theta_h)*N_h^(-theta_h) 
= (1-mu_h)*(1-N_h)^(-mu_h)*(C_h^mu_h*(1-N_h)^(1-mu_h))^(gamma_h-1);

/* intratemporal choice between C_f and N_f */
mu_f*C_f^(mu_f-1)*(C_f^mu_f*(1-N_f)^(1-mu_f))^(gamma_f-1)*omega_f*B_f^(sigma_f-1)*(omega_f*B_f^(sigma_f)+(1-omega_f)*A_f^(sigma_f))^((sigma_f^-1)-1)*exp(z_f)*K_f(-1)^theta_f*(1-theta_f)*N_f^(-theta_f) 
= (1-mu_f)*(1-N_f)^(-mu_f)*(C_f^mu_f*(1-N_f)^(1-mu_f))^(gamma_f-1);

/* intertemporal Eular equation on C_h */
mu_h*C_h^(mu_h-1)*(C_h^mu_h*(1-N_h)^(1-mu_h))^(gamma_h-1) 
= beta_h*mu_h*C_h(+1)^(mu_h-1)*(C_h(+1)^mu_h*(1-N_h(+1))^(1-mu_h))^(gamma_h-1)*(omega_h*A_h(+1)^(sigma_h-1)*(omega_h*A_h(+1)^(sigma_h)+(1-omega_h)*B_h(+1)^(sigma_h))^((sigma_h^-1)-1)*exp(z_h(+1))*theta_h*K_h^(theta_h-1)*N_h(+1)^(1-theta_h)+1-delta_h);

/* intertemporal eular equation on C_f */
mu_f*C_f^(mu_f-1)*(C_f^mu_f*(1-N_f)^(1-mu_f))^(gamma_f-1) 
= beta_f*mu_f*C_f(+1)^(mu_f-1)*(C_f(+1)^mu_f*(1-N_f(+1))^(1-mu_f))^(gamma_f-1)*(omega_f*B_f(+1)^(sigma_f-1)*(omega_f*B_f(+1)^(sigma_f)+(1-omega_f)*A_f(+1)^(sigma_f))^((sigma_f^-1)-1)*exp(z_f(+1))*theta_f*K_f^(theta_f-1)*N_f(+1)^(1-theta_f)+1-delta_f);

/* equal contribution of good A in h and f country */
mu_h*C_h^(mu_h-1)*(C_h^mu_h*(1-N_h)^(1-mu_h))^(gamma_h-1)*omega_h*A_h^(sigma_h-1)*(omega_h*A_h^(sigma_h)+(1-omega_h)*B_h^(sigma_h))^((sigma_h^-1)-1) 
= mu_f*C_f^(mu_f-1)*(C_f^mu_f*(1-N_f)^(1-mu_f))^(gamma_f-1)*(1-omega_f)*A_f^(sigma_f-1)*(omega_f*B_f^(sigma_f)+(1-omega_f)*A_f^(sigma_f))^((sigma_f^-1)-1);

/* equal contribution of good B in h and f country */
mu_h*C_h^(mu_h-1)*(C_h^mu_h*(1-N_h)^(1-mu_h))^(gamma_h-1)*(1-omega_h)*B_h^(sigma_h-1)*(omega_h*A_h^(sigma_h)+(1-omega_h)*B_h^(sigma_h))^((sigma_h^-1)-1) 
= mu_f*C_f^(mu_f-1)*(C_f^mu_f*(1-N_f)^(1-mu_f))^(gamma_f-1)*omega_f*B_f^(sigma_f-1)*(omega_f*B_f^(sigma_f)+(1-omega_f)*A_f^(sigma_f))^((sigma_f^-1)-1);

/* MC for good A */
A_h+A_f = exp(z_h)*K_h(-1)^theta_h*N_h^(1-theta_h);

/* MC for good B */
B_f+B_h = exp(z_f)*K_f(-1)^theta_f*N_f^(1-theta_f);

/* MC for final good in h country */
C_h+X_h = (omega_h*A_h^(sigma_h)+(1-omega_h)*B_h^(sigma_h))^(sigma_h^-1);

/* MC for final good in f country */
C_f+X_f = (omega_f*B_f^(sigma_f)+(1-omega_f)*A_f^(sigma_f))^(sigma_f^-1);

/* evolution of capital in h country */
X_h = K_h-(1-delta_h)*K_h(-1);

/* evolution of capital in f country */
X_f = K_f-(1-delta_f)*K_f(-1);

/* price of good A in terms of final good in h country */
qa_h = omega_h*A_h^(sigma_h-1)*(omega_h*A_h^(sigma_h)+(1-omega_h)*B_h^(sigma_h))^((sigma_h^-1)-1);

/* price of good B in terms of final good in h country */
qb_h = (1-omega_h)*B_h^(sigma_h-1)*(omega_h*A_h^(sigma_h)+(1-omega_h)*B_h^(sigma_h))^((sigma_h^-1)-1);

/* price of good A in terms of final good in f country */
qa_f = (1-omega_f)*A_f^(sigma_f-1)*(omega_f*B_f^(sigma_f)+(1-omega_f)*A_f^(sigma_f))^((sigma_f^-1)-1);

/* price of good B in terms of final good in f country */
qb_f = omega_f*B_f^(sigma_f-1)*(omega_f*B_f^(sigma_f)+(1-omega_f)*A_f^(sigma_f))^((sigma_f^-1)-1);

/* GDP of h country in terms of final good in h country */
Y_h = qa_h*exp(z_h)*K_h(-1)^theta_h*N_h^(1-theta_h);

/* GDP of f country in terms of final good in f country */
Y_f = qb_f*exp(z_f)*K_f(-1)^theta_f*N_f^(1-theta_f);

/* relative prices of foreign good in terms of home good - terms of trade */
p_h*qa_h = qb_h;
p_f*qb_f = qa_f;

/* net export - GDP ratio */
nx_h*Y_h = A_f-p_h*B_h;
nx_f*Y_f = B_h-p_f*A_f;

/* take log */
y_h = log(Y_h);
y_f = log(Y_f);
x_h = log(X_h);
x_f = log(X_f);
n_h = log(N_h);
n_f = log(N_f);
c_h = log(C_h);
c_f = log(C_f);
k_h = log(K_h);
k_f = log(K_f);
a_h = log(A_h);
a_f = log(A_f);
b_h = log(B_h);
b_f = log(B_f);

/* other terms */
ir_h = b_h-a_h;
ir_f = a_f-b_f;
rx_h*qa_f = qa_h;
rx_f*qb_h = qb_f;

/* TFP shocks */
z_h = rho_h_h*z_h(-1)+rho_h_f*z_f(-1)+eta_h;
z_f = rho_f_h*z_h(-1)+rho_f_f*z_f(-1)+eta_f;

end;
/* main model end */


@#for co in countries
beta_@{co} = 0.99;
mu_@{co} = 0.34;
gamma_@{co} = -1;
theta_@{co} = 0.36;
delta_@{co} = 0.025;
sigma_@{co} = -0.111;
@#endfor

omega_h = 0.873;
omega_f = 0.873;
rho_h_h = 0.97;
rho_h_f = 0.03;
rho_f_h = 0.03;
rho_f_f = 0.97;


initval;
@#for co in countries
C_@{co} = 0.4734;
N_@{co} = 0.5;
K_@{co} = 6.16;
X_@{co} = 0.154;
Y_@{co} = 0.6274;
nx_@{co} = 0;
ir_@{co} = -1.7347;
rx_@{co} = 1;
qa_@{co} = 0.5081;
qb_@{co} = 0.5081;
p_@{co} = 1;
z_@{co} = 0;
c_@{co} = log(0.4734);
n_@{co} = log(0.5);
k_@{co} = log(6.16);
x_@{co} = log(0.154);
y_@{co} = log(0.6274);

eta_@{co} = 0;
@#endfor

A_h = 1.0496;
A_f = 0.1852;
B_h = 0.1852;
B_f = 1.0496;
a_h = log(1.0496);
a_f = log(0.1852);
b_h = log(0.1852);
b_f = log(1.0496);
end;


shocks;
var eta_h; stderr 0.0073;
var eta_f; stderr 0.0073;
corr eta_h, eta_f = 0.290;
end;

steady;
check;

stoch_simul(order=1, hp_filter=1600 /*, nograph */) y_h c_h x_h n_h a_h b_h p_h nx_h rx_h    y_f c_f x_f n_f a_f b_f p_f nx_f rx_f;
