clear

num_threads_SL(6); % to use multiple cores

% Hamiltonian parameter
NC = 1;
T = 1e-12;
D = 1.;
J1 = 0.3;

% NRG parameter
Lambda = 3;
N = ceil(-6*log(T/10)/log(Lambda)); %% Lambda^(-N/2) ~ E_N ~ T / 10 ==> N ~ -2 log(T/10) / log(Lambda)
Nkeep = 1000;
% RhoV2in = [0.1;0.1]; % values outside of the 'ozin' grid are assumed to be zero
ozin = [-D; D];
RhoV2in = [1; 1];

% % % % U(1) charge * SU(2) spin
% Define operators
[FF,ZF,SF,IF] = getLocalSpace('FermionS','Z2charge','NC',NC);
% request a spinful-fermion(FermionS) local Hilbert space with abelian U(1) and non-Abelian SU(2) formed from NC channels
% FF: annihilator (spinor), ZF: fermionic parity (for anticommutation),
% SF: spin operator, f^dagger_alpha sigma_alpha,beta f_beta., IF: identity and other metadata


[Fs,Zs,Ss,Is] = setItag('s00','op',FF,ZF,SF,IF.E);
% operators Fs, Zs etc are now set to act on the `s00` bath zeroth site space.
[FL,ZL,SL,IL] = setItag('L00','op',FF,ZF,SF,IF.E);
% these operators instead are now set to act on the `L00` or impurity space. Only the spin SF 
% is retagged into SL, because the full spinor is not needed for the Kondo model.

%% Physical-spin operators: [S+/sqrt(2), S-/sqrt(2), Sz]

Sf = QSpace(1,3);
Sf(1) = quadOp(FL(1),FL(2),'*')/sqrt(2);
Sf(2) = quadOp(FL(2),FL(1),'*')/sqrt(2);
Sf(3) = makeIrop((quadOp(FL(1),FL(1),[]) ...
	- quadOp(FL(2),FL(2),[])))/2;
Sf = setItag('L00','op',Sf);

S0 = QSpace(1,3);
S0(1) = quadOp(Fs(1),Fs(2),'*')/sqrt(2);
S0(2) = quadOp(Fs(2),Fs(1),'*')/sqrt(2);
S0(3) = makeIrop((quadOp(Fs(1),Fs(1),[]) ...
	- quadOp(Fs(2),Fs(2),[])))/2;
S0 = setItag('s00','op',S0);

%% Charge-isospin operators: [C+/sqrt(2), C-/sqrt(2), Cz]

C0 = QSpace(1,3);

% C+ = c_up^dagger c_down^dagger
C0(1) = quadOp(Fs(1),Fs(2)','*')/sqrt(2);

% C- = c_down c_up
C0(2) = quadOp(Fs(2)',Fs(1),'*')/sqrt(2);

% Cz = (n_up + n_down - 1)/2
C0(3) = makeIrop((quadOp(Fs(1),Fs(1),[]) ...
	+ quadOp(Fs(2),Fs(2),[])...
	- Is))/2;

C0 = setItag('s00','op',C0);

iOdd = find(IL.Q{1}(:,1) == 1);
IL = getsub(IL,iOdd);
ZL = getsub(ZL,iOdd);

A0 = getIdentity(IL,2,Is,2,'K00*',[1 3 2]);

for J2 = [J1]..., 1.01 * J1,]
	nrgdata = sprintf('./_data/C2CK/%g-%g-%g-%d-%g', J1, J2, T, Nkeep, Lambda)
	% 2J ensures the term is S_f . sigma_0 and not S_0
	HSS = (2*J1) * (...
		contract(Sf(1),'!2*',S0(1), [2 1 3 4])...
		+ contract(Sf(2),'!2*',S0(2), [2 1 3 4])...
		+ contract(Sf(3),'!2*',S0(3), [2 1 3 4])...
	);

	% same significance of 2J
	HSC = (2*J2) * (...
		contract(Sf(1),'!2*',C0(1),[2 1 3 4])...
		+ contract(Sf(2),'!2*',C0(2),[2 1 3 4])...
		+ contract(Sf(3),'!2*',C0(3),[2 1 3 4])...
	);

	H0 = HSS + HSC;
	ops = [...
		contract(Fs(1),'!1',H0,[3 4 1 5 2]) - contract(H0,'!13',Fs(1));... Spec. func. for up
		contract(Fs(2),'!1',H0,[3 4 1 5 2]) - contract(H0,'!13',Fs(2));... spec. func. for down
		Sf(3); ... spin susceptibility
		]

	H0 = contract(A0,'!2*',{H0,'!13',A0}) + 1e-40*getIdentity(A0,2);

	ff = doZLD(ozin,RhoV2in,Lambda,N,1);

	NRG_SL(nrgdata,H0,A0,Lambda, ff{1}(2:end),FF,ZF, 'Nkeep',Nkeep);
	% getRhoFDM(nrgdata,T,'-v');
	%
	% [odisc,Adisc,sigmak] = getAdisc( ...
	% 	nrgdata, ...
	% 	ops(:),ops(:),ZF, ...
	% 	'Z_L00', ZL,...
	% 	'cflag', [1; 1; -1],...
	% 	'zflag', [1; 1; 0]...
	% 	);
	%
	% [ocont, Aup] = getAcont(odisc, Adisc{1}, sigmak, T/5);
	% [ocont, Adown] = getAcont(odisc, Adisc{2}, sigmak, T/5);
	% [ocont, ChiCont] = getAcont(odisc, Adisc{3}, log(Lambda), T/5);

	plotE( ...
		nrgdata, ...
		'title',sprintf( ...
			'$J_2 / J_1 = %g, T = %g, \\Lambda = %g, N_\\mathrm{keep} = %d$', ...
			J2 / J1,T,Lambda,Nkeep), ...
		'Emax',3);

	hflow = gcf;
	ax = gca;
	ax.Title.Interpreter = 'latex';
	ax.Title.FontSize = 16;
	exportgraphics( ...
		hflow, ...
		sprintf( ...
			'compact-2CK-NRG-J1=%g-J2=%g-Lambda=%g-N=%d-Nkeep=%d-T=%g.pdf', ...
			J1,J2,Lambda,N,Nkeep,T), ...
		'ContentType','vector', ...
		'BackgroundColor','white');

	close(hflow);

	% plotData = {0.5 * (Aup + Adown), ChiCont};
	% labels = {'$T(\omega)$', '$\chi^{\prime\prime}(\omega)$'};
	% fnames = {'SF', 'Chi'};
	% for op = 1:2
	% 	fig = figure('Color','white');
	% 	plot(ocont,plotData{op},'LineWidth',1.8);
	%
	% 	set(gca,'XScale','log');
	% 	if op == 2
	% 		set(gca,'YScale','log');
	% 	end
	%
	% 	ax = gca;
	% 	title(ax, ...
	% 		sprintf( ...
	% 			'$J_1=%g,\\quad J_2=%g, \\quad T=%g, \\quad \\Lambda=%g,\\quad N_{\\mathrm{keep}}=%d$', ...
	% 			J1,J2,T,Lambda,Nkeep), ...
	% 		'Interpreter','latex', ...
	% 		'FontSize',16);
	% 	grid on;
	%
	% 	xlabel('$\omega$','Interpreter','latex');
	% 	ylabel(labels{op},'Interpreter','latex');
	%
	% 	exportgraphics( ...
	% 		fig, ...
	% 		sprintf( ...
	% 			'compact-2CK-%s-J1=%g-J2=%g-T=%g-Lambda=%g-Nkeep=%d.pdf', ...
	% 			fnames{op},J1,J2,T,Lambda,Nkeep), ...
	% 		'ContentType','vector', ...
	% 		'BackgroundColor','white');
	%
	% 	close(fig);
	% end
end
