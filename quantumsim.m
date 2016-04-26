classdef quantumsim < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %constants, lenth in nm, energy in eV
        % # of grid points
        a=0.5;
        E_f=0.15;       % fermi energy in eV
        E_g=1;        % band gap in eV
        V_ds=0.5;       % drain-source voltage in V
        V_g=2;          % gate potential Psi_g=-e*V_g in eV
        d_ox=5;         % oxide thickness in nm
        d_ch=5;        % channel thickness in nm
        e=util.const.e; % elementary charge
        k_0=8.85e-12;   % dielectric constant
        k_Si=11.2;     % dielectric constant for Silicon
        k_ox=3.9;       % dielectric constant oxide
        geo=1;          % Geometriefaktor # of gates, 'w' wrapgate
        l_ch=40;        % channel length
        l_ds=40;        % length of drain and source regions
        N_dot=0;        % # Dopands
        lambda
        L_sparse
        n_left
        n_right
        N
        Psi_g
        Psi_bi
        rho
        Psi_f
        Psi_0
        epsilon=10e-15;         % tolarance for fermi function
        E_fs=0.05;  %?????
        E_fd
        T=300;  %temperature in Kelvin
        dE=0.001; % energy step
        
        E_max
        E
        % From here: (Eigenenergies)
        H
        G_r
    end
    
    properties (Dependent)
        
    end
    
    methods
        function self =  quantumsim() %constructor
            self.get_lambda();
            self.get_N();
            self.get_laplacian();
            self.E_fd=-self.V_ds+0.05;
            self.E_max=self.E_fs-util.const.k_b*self.T*log(self.epsilon)/util.const.e;
            self.E=self.Psi_0:self.dE:self.E_max;
        end
        
        function get_lambda(self)
            self.lambda=sqrt(self.k_Si/self.k_ox*self.d_ch*self.d_ox/self.geo);
        end
        
        function get_laplacian(self)
            %tic;
            super=zeros(self.N,1);
            super(2:end)=1;
            super(2)=2;
            sub=zeros(self.N,1);
            sub(1:end-1)=1;
            sub(end-1)=2;
            middle=zeros(self.N,1);
            middle(:)=-2;
            self.L_sparse=spdiags([super,middle,sub],[1,0,-1],self.N,self.N)/self.a.^2;
            %toc;
        end
        
        function get_N(self)
            l_g=2*self.l_ds+self.l_ch;
            self.N=floor(l_g/(self.a))+1;
            
            self.n_left=floor(self.l_ds/self.a);      % source region 1 to n-left
            self.n_right = self.N - self.n_left;       % drain region n_right-end
            
        end
        
        function init_vectors(self)
            %create empty arrays
            self.Psi_g=zeros(self.N,1);
            self.Psi_bi=zeros(self.N,1);
            
            % soure both 0
            
            % channel
            self.Psi_g(self.n_left+1:self.n_right-1) =-self.V_g;
            self.Psi_bi(self.n_left+1:self.n_right-1)=self.E_f+self.E_g/2.;
            
            %drain
            self.Psi_g(self.n_right:end) =0;
            self.Psi_bi(self.n_right:end)=-self.V_ds;
            
            % create charge density
            self.rho=zeros(self.N,1);
        end
        
        function calc_potential(self)
            %tic;
            self.Psi_f=(self.L_sparse-spdiags(zeros(self.N,1)+1/self.lambda^2,0,self.N,self.N))\...
                ((self.rho+self.N_dot)/self.k_0/self.k_Si-1/self.lambda^2*(self.Psi_g+self.Psi_bi));
            %toc;
            self.Psi_0=max(self.Psi_f);
            self.E=self.Psi_0:self.dE:self.E_max;
        end
        
        function I = calc_current(self)
            f_s = @(E) 1./(exp((E-self.E_fs).*util.const.e./util.const.k_b./self.T)+1);
            f_d = @(E) 1./(exp((E-self.E_fd).*util.const.e./util.const.k_b./self.T)+1);
            I=2*util.const.e/util.const.h*sum(f_s(self.E)-f_d(self.E))*self.dE*util.const.e*1e-3;
        end
        
        function plot_potential(self)
            figure, plot((0:self.N-1).*self.a,self.Psi_f);
        end
        
        function set_V_ds(self,V)
            self.V_ds=V;
            self.E_fd=-self.V_ds+0.05;
            self.init_vectors();
        end
        
        function set_V_g(self,V)
            self.V_g=V;
            self.init_vectors();
        end
        
        function set_l_ch(self,l)
            self.l_ch=l;
            self.get_N();
            self.get_laplacian();
            self.init_vectors();
            
        end
        
        function S = plot_Vg_I(self,min,max,step,plotme)
            I=[];
            for V=min:step:max
                self.set_V_g(V);
                self.calc_potential();
                I(end+1)=self.calc_current();
            end
            x=min:step:max;
            y=log10(I);
            index = (x >= 0) & (x <= 0.4);
            p = polyfit(x(index),y(index),1);  %# Fit polynomial coefficients for line
            yfit = p(2)+x.*p(1);                %# Compute the best-fit line
            
            if nargin > 4
            figure, plot(x,y);                      %# Plot the data
            hold on;              %# Add to the plot
            plot(x,yfit,'r');     %# Plot the best-fit line           
            end
            S=1/p(1);
        end
        
        function plot_Vds_I(self,min,max,step)
            I=[];
            for V=min:step:max
                self.set_V_ds(V);
                self.calc_potential();
                I(end+1)=self.calc_current();
            end
            x=min:step:max;
            y=I;                    
            figure, plot(x,y);                                      
        end
        
        function plot_S(self)
            S=[];
            for l=15:5:100
                set_l_ch(self,l);
                S(end+1)=self.plot_Vg_I(0,1,0.01);
            end
            figure,plot(15:5:100,S);
        end
        
        %task 3
        function task3(self)
            
            self.set_l_ch(40);
            %check other stuff!
            Vg_min=0;
            Vg_max=1;
            Vg_step=0.1;
            
            Vd_min=0;
            Vd_max=1;
            Vd_step=0.01;
            
            figure();
            hold on;
            
            for Vd=Vg_min:Vg_step:Vg_max
            self.set_V_g(Vd);
            I=[];
            for V=Vd_min:Vd_step:Vd_max
                self.set_V_ds(V);
                self.calc_potential();
                I(end+1)=self.calc_current();
            end
            x=Vd_min:Vd_step:Vd_max;
            y=I;                          
            plot(x,y);                      %# Plot the data
                         %# Add to the plot                       
            end
            hold off;
            
            self.plot_S();
            
        end
        
        function calc_green(self)
            self.E=0:0.000001:0.0004;
            super=zeros(self.N,1);
            super(2:end)=1;
            sub=zeros(self.N,1);
            sub(1:end-1)=1;
            middle=zeros(self.N,1);
            middle(:)=-2;
            pot = zeros(self.N,1);
            eta = 0.0000008*(1i);
            m=0.9*util.const.m_e;
            t = util.const.h_bar^2/(2*m*(self.a*10^(-9))^2*util.const.e);
            self.H=-t.*spdiags([super,middle,sub],[1,0,-1],self.N,self.N)+spdiags(pot,0,self.N,self.N);
            self.G_r= zeros(length(self.E),self.N);
            for k=1:length(self.E)
                self.G_r(k,:) = imag(diag(inv(spdiags(ones(self.N,1)*...
                    self.E(k)+eta,0,self.N,self.N)-self.H)/self.a)); 
            end
            figure,imagesc(self.G_r);
            colorbar();
            set(gca,'Ydir','Normal');
            
            disp(-1/pi*sum(self.G_r(:))*self.a*0.000001)
        end
        
        
        
        
    end
    
end

