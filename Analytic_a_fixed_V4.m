%% Analytic_a_fixed.m
% MPAdS: Unified Angular Eigenvalues + Analytic QNM only
% Scan over ell while fixing a (a_fixed). Calculates A_jlm, Cho coeffs f1,f2,
% analytic QNM omega and a physically consistent analytic_status including
% superradiant vs non-superradiant instabilities.

clear; clc;
fprintf('=== MPAdS: Unified Angular Eigenvalues + Analytic QNM only ===\n');

%% ------------------ User Parameters ------------------
D = 4; n = D-4;

 radial_n = 1;j = 0;l = 3; m = 1;

a_fixed = 0.3;                % fixed rotation parameter a
ell_list = 1:0.2:700;       %radial_n = 0;j = 0;l = 3; m = 3;


mu_s = 0.045;
N_cheb = 300;
save_prefix = 'MPAdS_AnalyticOnly_EllScan';
               % radial mode index

G = 1;              % Newton constant
M = 1;        % ---->Blackhole MASS

a_cutoff = 1000;             % cutoff used for Cho analytic validity (if used)

%% ------------------ Preallocate ------------------
A_vals = zeros(size(ell_list));
F1_vals = NaN(size(ell_list));
F2_vals = NaN(size(ell_list));
omega_vals = NaN(size(ell_list));
status_vals = strings(size(ell_list));
rplus_vals  = NaN(size(ell_list));
rminus_vals = NaN(size(ell_list));
rm_vals = NaN(size(ell_list));
%% ------------------ Compute <cos^2> expectation ------------------
cos2_expect = compute_cos2_expectation(D,l,N_cheb);
SD2 = 2*pi^((D-1)/2) / gamma((D-1)/2);
%% ------------------ Main Loop (scan ell) ------------------


for k = 1:length(ell_list)
      % ---------- Physical constraint on j ----------
if D == 4 && j ~= 0
    fprintf('Skipping j = %d because j > 0 is not physical in 4D.\n', j);
    continue;   % j > 0 not physical in 4D
end


             % Condition 1: l >= j + Abs[m]
          if l < (j + abs(m)) 
              fprintf('Skipping m = %d because l < j + Abs[m] (l = %d, j = %d)\n', m, l, j);
            continue; % Skip the rest of the loop for this m
          end
         if radial_n ~= 0.5 * (l - (j + abs(m)))
          fprintf('Skipping j = %d, m = %d because (radial_n = %d,l = %d)\n', j, m, radial_n,l);
          continue;
         end

                %fprintf('l=%d, m=%d, j=%d, n_r=%d\n', l, m, j, radial_n);
    ell = ell_list(k);
    a_val = a_fixed;
omega_n = (1/ell) * (2*radial_n + l + (n+3)/2 + sqrt(((n+3)/2)^2 + mu_s^2 * ell^2));
omega=omega_n;
tilde_omega = sqrt(max(0, omega^2 - mu_s^2));
     alpha = - (a_val^2 / ell^2);      % Cho/AdS convention
    c = a_val * tilde_omega;
    %----- Sigma(a)
    Sigma = 1 - (a_val^2)/(ell^2);

      %----- Compute rm from fixed M = 1
    numer = 8*pi*G * Sigma^2*M;
    denom = SD2 * (1 + (D-4)/2 * Sigma);
    r_m = ( numer / denom )^(1/(D-3));   % <--- NEW
   rm_vals(k)=r_m;
% --- 2) Radial horizons ---
    [r_plus, r_minus] = compute_horizons(a_val, ell, r_m, n);
    rplus_vals(k)  = r_plus;
    rminus_vals(k) = r_minus;
 % --- 1) Angular eigenvalue ---
    %if m==0 && a_val <= a_cutoff
     if a_val <= a_cutoff
        [f1,f2,A0] = f_coeffs_Cho_AdS_analytic(D,l,m,j,alpha,c);
        delta_f2_mass = 0.5 * mu_s^2 * cos2_expect;
        f2 = f2 + delta_f2_mass;
        A_jlm = A0 +f1*a_val*tilde_omega + f2*(a_val*tilde_omega)^2;
        F1_vals(k) = f1;
        F2_vals(k) = f2;
        else
    if m==0 
        A_jlm=NaN;
        F1_vals(k) = NaN;
        F2_vals(k) = NaN;
    else
        A_jlm = compute_A_jlm(a_val,N_cheb,tilde_omega,mu_s,ell,D,l,m,j);
        F1_vals(k) = NaN;
        F2_vals(k) = NaN;
    end
      
eps=A_jlm+(mu_s^2*r_plus^2)-l*(l+n+1);
     end
    A_vals(k) = A_jlm;
  
   
    Omega_plus = a_val*(1 - a_val^2/ell^2)/(r_plus^2 + a_val^2);

    F1 = r_plus^(n+1) - r_minus^(n+1);
          cosh_term3 = cosh(2*pi * r_plus^(2+n) * (omega -(m* Omega_plus)) / ((1+n) * F1));
    sinh_term3 = sinh(2*pi * r_plus^(2+n) * (omega -(m*Omega_plus)) / ((1+n) * F1));
    gamma_arg3 = sym(1+l/(1+n) - 2i * r_plus^(2+n) * (omega - (m*Omega_plus)) / ((1+n) * F1));
  abs_gamma3 = double(abs(gamma(gamma_arg3))^2);
       sqrt_term = sqrt(((3+n)/2)^2 + ell^2 * mu_s^2);

  eta= (( F1^(1 + (2*l)/(1 + n)) * abs_gamma3 * ...
          gamma(1 + l/(1 + n))^2 * gamma((3 + n)/2 +l+radial_n+ sqrt_term) * prod( l + (n+1)/2 ...
          + (1:radial_n) )) / ( pi * ell^(2 *l + n + 2)* gamma(radial_n + 1)*gamma(1 + (2*l)/(1+n)) * gamma(2 + (2*l)/(1+n)) ...
        * gamma(l + (1+n)/2) * gamma(1 +radial_n+ sqrt_term)));
       
        re_omega = (2*radial_n + (3 + n)/2+l + sqrt_term)/ell - eta *(cosh_term3* tan(pi/(1+n) * (l + eps/(2*l + n + 1))));

        im_omega =  - eta * sinh_term3; 

        omega_new = re_omega + 1i*im_omega;
      
  
    omega_vals(k) = omega_new;

 % --- Determine analytic_status (physically consistent) ---
    if imag(omega_new) < 0
        analytic_status = "STABLE";

    elseif imag(omega_new) > 0
        % unstable: decide if superradiant
        if m > 0 && real(omega_new) > 0 && real(omega_new) < m * Omega_plus
            analytic_status = "SUPERRADIANT"; % superradiant instability
         elseif m < 0 && real(omega_new) < 0 && real(omega_new) < m * Omega_plus
            analytic_status = "SUPERRADIANT"; % superradiant instability
        else
            analytic_status = "NON_SUPERRADIANT_UNSTABLE"; % other instabilities
        end

    else
        analytic_status = "MARGINAL"; % Im == 0
    end
    status_vals(k) = analytic_status;
end


%% ------------------ Save CSV ------------------
T0 = table((ell_list)',(a_fixed ./ ell_list)', A_vals', F1_vals', F2_vals', ...
    real(omega_vals)', imag(omega_vals)', ...
    rplus_vals', rminus_vals',rm_vals', status_vals', ...
    'VariableNames', {'ell','a/ell','A_jlm','f1','f2', ...
    'omega_real','omega_imag','r_plus','r_minus','rm','analytic_status'});

%writetable(T,[save_prefix '_results_a_fixed.csv']);
%fprintf('Saved results to %s_results_a_fixed.csv\n',save_prefix);
writetable(T0,[save_prefix '_results_a_fixed_mu0.csv']);
disp(T0);

%% ================== Local Functions ==================

function [r_plus, r_minus] = compute_horizons(a_val, ell, r_m, n)
    % Define Delta(r)
    Delta = @(r) (r.^2 + a_val.^2) .* (1 + r.^2 ./ ell.^2) ...
                 - r_m^(n+1) ./ r.^(n-1);

    % Find roots using fzero on interval grid
    r_grid = linspace(1e-6, 50*ell, 15000);
    Dvals = Delta(r_grid);

    % Locate sign changes
    idx = find(diff(sign(Dvals)));

    roots_list = [];
    for k = 1:length(idx)
        r1 = r_grid(idx(k));
        r2 = r_grid(idx(k)+1);
    
       try 
       rroot = fzero(Delta, [r1, r2]);
           if rroot > 0
                roots_list(end+1) = rroot; %#ok<AGROW>
           end
       catch ME
       roots_list = NaN;
               fprintf('ERROR at roots=%.3f: %s\n',roots_list,ME.message);

       end
    end

    roots_list = sort(unique(round(roots_list, 10)));

    if isempty(roots_list)
        r_plus = NaN;
        r_minus = NaN;
        return;
    end

    % Outer and inner horizons
    if isscalar(roots_list)
        r_plus = roots_list(1);
        r_minus = 0;
    else
        r_minus = roots_list(1);
        r_plus  = roots_list(end);
    end
end
    %r_plus = max(roots_list);
   % if length(roots_list) >= 2
       % r_minus = min(roots_list);
   % else
      %  r_minus = 0;   % fallback
   % end
%end

function val = compute_cos2_expectation(D,l,N)
    n = D-4;
    [x,~] = cheb(N);
    theta = (x+1)*(pi/4);
    w_full = sin(theta).*cos(theta).^n;
    w_full(abs(w_full)<1e-14) = 1e-14;
    Y0 = sqrt(w_full) .* legendreP(l, cos(theta));
    normY0 = sqrt(trapz(theta, Y0.^2 .* w_full));
    if normY0==0, normY0=1; end
    Y0 = Y0 / normY0;
    val = trapz(theta, Y0.^2 .* cos(theta).^2 .* w_full);
end

function [f1,f2,A0] = f_coeffs_Cho_AdS_analytic(D,l,m,j,alpha,c)
    n = D-4; A0 = l*(l+n+1); m2 = m^2;
    P1 = 2*l^4 + 4*(n+1)*l^3 + (2*m2 + 3*n*(n+2) -1)*l^2 + (n+1)*(2*m2+n^2+2*n-3)*l + m2*(n+1)*(n+3)...
    - j^2*(n^2+2*l*n+4*n+2*l*(l+1)+3) - j*(n-1)*(n^2+2*l*n+4*n+2*l*(l+1)+3);
    den_f1 = (2*l+n-1)*(2*l+n+3); if abs(den_f1)<1e-14, den_f1=den_f1+1e-14; end
    if abs(c)<1e-14, f1=0; else 
        f1 = (alpha/c)*(P1/den_f1);
    end
    % compute f2 using original formula with small-c safeguard
     A1 = 2*l^2 + 2*l*(n-1) + n*(n+2) + 5;
    A2 = 2*l^2 + 2*l*(n+3) + (n+3)^2;
    t1 = (j-l-abs(m))*(j-l+abs(m))*(j+l-abs(m)+n-1)*(j+l+abs(m)+n-1);
    t2 = (j-l-abs(m)-2)*(j-l+abs(m)-2)*(j+l-abs(m)+n+1)*(j+l+abs(m)+n+1);
    den1 = (2*l+n-3)*(2*l+n-1)^3*(2*l+n+1);
    den2 = (2*l+n+1)*(2*l+n+3)^3*(2*l+n+5);
    if abs(den1)<1e-16, den1=den1+1e-16; end
    if abs(den2)<1e-16, den2=den2+1e-16; end
    alpha_part = (-c^2*A1*t1/den1) + (c^2*A2*t2/den2);
    PolyBig = -2*j^2*(l^2+(l-1)*n+l-3) - 2*j*l*(n+3)*(l+n+1) - 2*l^4 -4*l^3*(n+1) + l^2*(2*m2-3*n*(n+2)+1) - l*(n+1)*(-2*m2+n^2+2*n-3)-2*m2*(n+3);
    pref_poly = 2/((2*l+n-1)*(2*l+n+3));
    term_alpha_part = alpha_part + pref_poly*PolyBig + 2*j*(j+n+1) - 2*m2;
    alpha_bracket = 0.5*alpha*term_alpha_part;
    inner_frac = (-j^2 - j*n + j + l*(l+n+1) + m2 + n -1)/((2*l+n-1)*(2*l+n+3));
    c_term = c*(f1 - 2*c*inner_frac);
    big_bracket = alpha_bracket + c_term;
    
        f2 = - (1 / c^2) * big_bracket;
    end

function A_jlm = compute_A_jlm(a_val,N,omega,mu_s,ell,D,l,m,j)
    n = D-4;
    [x,Dcheb] = cheb(N);
    theta = (x+1)*(pi/4);
    dth_dx = pi/4; D_theta = Dcheb/dth_dx;
    Delta_theta = 1 - (a_val^2)/(ell^2)*cos(theta).^2;
    Sigma = 1 - (a_val^2)/(ell^2);
    w_full = sin(theta).*cos(theta).^n; w_full(abs(w_full)<1e-14)=1e-14;
    sin_th = sin(theta); sin_th(abs(sin_th)<1e-12)=1e-12;
    cos_th = cos(theta); cos_th(abs(cos_th)<1e-12)=1e-12;
    Pot_full = - j*(j+n-1)./cos_th.^2 - (1./Delta_theta).*((m*Sigma./sin_th - a_val*omega.*sin_th).^2) + (a_val^2)*(mu_s^2).*cos_th.^2;
    ind = 2:N; W = w_full(ind); Delta_int = Delta_theta(ind); Pot = Pot_full(ind); Dth_int = D_theta(ind,ind);
    A_op = Dth_int*(diag(Delta_int.*W)*Dth_int);
    L_int = diag(1./W)*A_op + diag(Pot);
    lam = eig(L_int); A_sph = l*(l+n+1); [~,idx] = min(abs(lam + A_sph)); A_jlm = -lam(idx);
end

function [x,D] = cheb(N)
    if N==0, x=1; D=0; return; end
    k = (0:N)'; x = cos(pi*k/N); c = [2; ones(N-1,1);2].*(-1).^k;
    X = repmat(x,1,N+1); dX = X - X'; D = (c*(1./c)')./(dX+eye(N+1)); D = D - diag(sum(D,2));
end







   
%% ============================================================
%  Im(omega) vs a/ell with vertical region boundaries
% ============================================================
%a_fixed./
% ---------- Data ----------
c1 = 3e8;           % speed of light [m/s]
G1 = 6.674e-11;     % Newton constant [m^3/kg/s^2]
M1_sun = 1.989e30;  % solar mass [kg]
r_g=G1*M1_sun/c1^2;
 Mphys = 10 * M1_sun;
omega_SI_real = real(omega_vals) * (c1^3 / (G1 * Mphys));   % rad/s
omega_SI_imag = imag(omega_vals) * (c1^3 / (G1 * Mphys));

y = 1./abs(omega_SI_imag / (2*pi))/(60*60*24*365);   % Hz

x = 10^-3*10*r_g*ell_list;        % a/ell

% Remove NaNs (safety)
valid = ~(isnan(x) | isnan(y));
x = x(valid);
y = y(valid);
status = status_vals(valid);

%% ---------- Figure ----------
figure; hold on; box on;


% Main curve
plot(x, y,'k-', 'LineWidth', 1);
set(gca,'yDir','reverse');
%set(gca,'xDir','reverse');
xlabel('$\ell (Km)$','Interpreter','latex');
%ylabel('$\tau (years)$','Interpreter','latex');
ylabel('$\tau\,\, (years)$','Interpreter','latex');

title(['($jnlm$ = ',num2str(j),num2str(radial_n),num2str(l),num2str(m),')'], ...
      'Interpreter','latex');

grid off;



%% ---------- Vertical boundaries: Superradiant region ----------
idx_sr = status == "SUPERRADIANT";

if any(idx_sr)
    x_sr_min = min(x(idx_sr));
    x_sr_max = max(x(idx_sr));

    xline(x_sr_min,'g','LineWidth',1,...
          'Label','SR onset',...
          'Interpreter','latex');

   % xline(x_sr_max,'g','LineWidth',1,...
      %'Label','SR end',...
        %  'Interpreter','latex');
end

%% ---------- Point classification (visual aid) ----------
idx_stable = status == "STABLE";
idx_nsr    = status == "MARGINAL";

plot(x(idx_stable), y(idx_stable), ...
     'ro','MarkerFaceColor','g','MarkerSize', 2, 'LineWidth', 1);

plot(x(idx_sr), y(idx_sr), ...
     'bo','MarkerFaceColor','r','MarkerSize', 2, 'LineWidth', 1);

plot(x(idx_nsr), y(idx_nsr), ...
     'ko','MarkerFaceColor','k','MarkerSize', 2, 'LineWidth', 1);

%% ---------- Legend ----------
legend({'$\tau=\frac{1}{Im(\omega)}$', ...
        'SR onset', ...
        'Stable','Superradiant','Non-SR unstable'}, ...
       'Interpreter','latex','Location','best');



ylim([0,0.7*10^8]);%radial_n = 1;j = 0;l = 3; m = 1;
