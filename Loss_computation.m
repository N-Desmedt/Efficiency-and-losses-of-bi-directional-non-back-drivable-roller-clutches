clear; close all
mu=0.244; mu_k=0.31;nu_f=0.44; nu_fk=0.44; alpha=6.5*pi/180;m=0.032;R=0.0255;r=R/3;F_Nk=2.3;N=4;
c2a=cos(2*alpha); s2a=sin(2*alpha); %precomputation for simplification
for i=1:2500
    omega=i/60*2*pi; %speed range definition    
    F_c=m*(R-r)*omega^2; %Roller entrifugal force (equal for all state)
    %% state 1
    A=-mu*mu_k*F_Nk*((mu+1)*s2a+(mu-1)*c2a+mu);
    B=F_Nk*((mu_k-mu^2)*s2a+(mu*mu_k+mu)*c2a)-mu*F_c*(mu+mu*c2a)-mu*R*omega*nu_f*(1+c2a);
    C=-mu*F_Nk-mu*F_c*s2a-R*omega*nu_f*s2a;
    betap=(-B+sqrt(B^2-4*A*C))/2/A; betam=(-B-sqrt(B^2-4*A*C))/2/A; % compute the two possible value of the parameter beta
    if (betap>0 && betap<1)||(betam>0 && betam<1) %check condition state 1 
        if (betap>0 && betap<1) beta=betap; else beta=betam; end %select right solution for beta
        F_fk=mu_k*beta*F_Nk;
        F_Nr=(F_Nk+F_c*s2a-R*omega*nu_f*(c2a+1)+F_fk)/(s2a+mu*c2a+mu);
        F_fr=mu*F_Nr+nu_f*R*omega;
        P_loss(i)=N/2*(F_fr*R*omega); % compute losses 
    end
    %% state 2
    omega_r=((mu^2*mu_k+mu+(mu_k*mu^2-mu-2*mu_k*mu)*c2a+(mu_k*mu^2+mu^2+mu*mu_k-mu_k)*s2a)*F_Nk+(mu^2*c2a+mu*s2a+mu^2)*F_c+(s2a+mu*c2a+mu)*nu_f*R*omega)/(((2*s2a+2*mu+2*mu*c2a)*nu_f-((mu^2+mu-1)*s2a+(mu^2-2*mu)*c2a+mu^2)*nu_fk)*r);
    v_r=R*omega-r*omega_r; v_k=r*omega_r; v_s=v_k;
    F_Ns=(((1-mu_k*mu)*c2a-(mu+mu_k)*s2a)*F_Nk-mu*F_c-nu_f*v_r-(s2a+mu*c2a)*nu_fk*v_k-(c2a-mu*s2a)*nu_f*v_s)/(s2a+2*mu*c2a-mu^2*s2a);
    if omega_r<omega*R/r && F_Ns>0 %check condition state 2
        F_Nr=((mu*mu_k+1)*F_Nk+(mu*c2a+s2a)*F_c-(c2a-mu*s2a)*nu_f*v_r+mu*nu_fk*v_k-nu_f*v_s)/(s2a+2*mu*c2a-mu^2*s2a);
        F_fs=mu*F_Ns+nu_f*v_s; F_fr=mu*F_Nr+nu_f*v_r; F_fk=mu_k*F_Nk+nu_fk*v_k;
        P_loss(i)=N/2*(F_fs*v_s+F_fr*v_r+F_fk*v_k); % compute losses 
    end
    %% state 3:
    F_Ns=((c2a-mu_k*s2a-mu_k)*F_Nk-(nu_f+nu_f*c2a+nu_fk*s2a+nu_fk)*R*omega)/((s2a+mu+mu*c2a));
    F_Nr=((1-mu_k*c2a)*F_Nk-mu*(1+c2a)*F_Ns+F_c*s2a-(nu_f*c2a+nu_f+nu_fk*c2a)*R*omega)/s2a;
    F_fr=mu*F_Ns+mu_k*F_Nk+(nu_f+nu_fk)*R*omega;
    if F_fr<mu*F_Nr&& F_Ns>0 %check condition state 3
        v_s=R*omega; v_k=v_s;
        F_fs=mu*F_Ns+nu_f*v_s; F_fk=mu_k*F_Nk+nu_fk*v_k;
        P_loss(i)=N/2*(F_fs*v_s+F_fk*v_k); % compute losses 
    end
    %% state 4
    omega_r=((mu_k*mu*s2a+mu*mu_k-mu*c2a)*F_c-(c2a-mu_k-mu_k*s2a)*nu_f*R*omega)/(((mu-c2a+mu*s2a)*nu_fk+(mu_k-c2a+mu_k*s2a)*nu_f)*r);
    F_Nr=((nu_fk+nu_f-mu_k*nu_f*c2a)*r*omega_r-(1-mu_k*c2a)*nu_f*R*omega-mu_k*F_c*s2a)/((mu-mu_k*s2a-mu*mu_k*c2a));
    F_Nki=(mu*F_Nr-(nu_fk+nu_f)*r*omega_r+nu_f*R*omega)/mu_k;
    if omega_r<omega*R/r && F_Nki>F_Nk %check condition state 4
        v_r=R*omega-r*omega_r; v_k=r*omega_r;
        F_fk=mu_k*F_Nki+nu_fk*v_k; F_fr=mu*F_Nr+nu_f*v_r;
        P_loss(i)=N/2*(F_fr*v_r+F_fk*v_k); % compute losses 
    end
    %% 2-4 transition state
    omega_r=((mu*mu_k+(mu_k+mu+mu*mu_k)*s2a+(mu*mu_k-1-mu_k)*c2a)*	F_Nk+(mu*c2a+mu)*F_c+(1+c2a)*nu_f*R*omega)/(-((mu+(mu+1)*s2a+(mu-1)*c2a)*nu_fk-(1+c2a)*nu_f)*r);
    F_Nr=(	(1+mu_k)*F_Nk-(nu_f+nu_f*c2a)*R*omega+(nu_fk+nu_f+nu_f*c2a)*r*omega_r+F_c*s2a)/((mu+s2a+mu*c2a));
    gamma=(mu*F_Nr+nu_f*R*omega-mu_k*F_Nk-(nu_fk+nu_f)*r*omega_r)/(nu_f*r*omega_r);
    if omega_r<omega*R/r && gamma>0 && gamma <1 %check condition transition state 2-4
        v_s=r*omega_r; v_r=R*omega-r*omega_r; v_k=v_s;
        F_fs=gamma*nu_f*v_s; F_fk=mu_k*F_Nk+nu_fk*v_k; F_fr=mu*F_Nr+nu_f*v_r;
        P_loss(i)=N/2*(F_fs*v_s+F_fr*v_r+F_fk*v_k); % compute losses 
    end
    %% state 5
    omega_r=omega*R/r; v_k=r*omega_r;
    F_Nki=((s2a+1)*nu_fk*v_k)/(-mu_k-mu_k*s2a+c2a);
    F_Nr=(F_c*s2a-nu_fk*v_k*c2a+(1-mu_k*c2a)*F_Nki)/s2a;
    F_fr=mu_k*F_Nki+nu_fk*v_k;F_fk=mu_k*F_Nki+nu_fk*v_k;
    if F_fr<mu*F_Nr && F_Nki>F_Nk %check condition state 5
        P_loss(i)=N/2*(F_fk*v_k); % compute losses 
    end
    %% 3-5 transition state
    gamma=(-(mu_k*s2a+mu_k-c2a)*F_Nk-(s2a+1)*nu_fk*R*omega)/((c2a+1)*nu_f*R*omega);
    F_Nr=((1-mu_k*c2a)*F_Nk-(1+c2a)*gamma*nu_f*R*omega+F_c*s2a-nu_fk*R*omega*c2a)/s2a;
    F_fr=mu_k*F_Nk+nu_fk*R*omega+gamma*nu_f*R*omega;
    if F_fr<mu*F_Nr && gamma>0 && gamma <1 %check condition transition state 3-5
        v_k=R*omega; v_s=v_k;
        F_fs=gamma*nu_f*v_s; F_fk=mu_k*F_Nk+nu_fk*v_k;
        P_loss(i)=N/2*(F_fk*v_k+F_fs*v_s); % compute losses 
    end
end
