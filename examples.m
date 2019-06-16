%% Calculate Wigner function of a single mode squeezed vacuum

C = create_vacuum(1);
C = squeeze(C,1); %Squeezing factor is 1
W = calculate_single_mode_wigner(C,100,0.1); %Wigner function calculated on a 100x100 grid with step 0.1
surf(W)

%% Calculate photon number statistics of a coherent state with average photon number of 1

C = create_vacuum(1);
C = displace(C,sqrt(2)*1); %displacement = 1 in the x and p basis, which is equivalent to alpha=sqrt(2)
P = photon_number_stats(C,10); %Calculate statistics up to 10 photons
bar3(P)

%%
% Calculate photon number statistics of a 2 mode squeezed vacuum, created by 
% interfering two single mode squeezed vacuums on a 50:50 beam splitter 

C = create_vacuum(2);
C = squeeze(C,[1,1]); %Squeezing factor for each mode is 1
C = beam_splitter(C,[1,2],pi/4); %beam splitter op between modes 1 and 2, with reflectivity cos(pi/4)
P = photon_number_stats(C,10); %Calculate statistics up to 10 photons in each mode
bar3(P)

%% Calculate photon number statistics of a complex two mode state, to showcase all ops

C = create_vacuum(2);
C = amplify(C,[0,0.7]); %Phase-insensitive amplification on mode 2
C = squeeze(C,[1,0.5]);
C = phase_shift(C,[pi/8,-pi/3]); 
C = displace(C,[0.5+2*1i,0]); %Displace mode 1
C = beam_splitter(C,[1,2],pi/3); %beam splitter op between modes 1 and 2, with reflectivity pi/3
C = add_loss(C,[0,0.2]); %loss added to mode 2
P = photon_number_stats(C,10); %Calculate statistics up to 10 photons in each mode
bar3(P)

