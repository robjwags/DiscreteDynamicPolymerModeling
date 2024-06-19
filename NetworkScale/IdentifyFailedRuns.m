function package = IdentifyFailedRuns(BSOM,Controls)

global length_conversion damper_conversion...
    Np Nt Ns b D dt...
    samples eq_time_factors N_Kuhns phis kds Weissenbergs...
    sample eq_time_factor N_Kuhn phi kd Weissenberg...
    output_folder OutputAtom_rlx OutputAtom...
    BeadSpringOrMeso LinearOrLangevin

BeadSpringOrMeso = BSOM;
LinearOrLangevin = Controls.LinearOrLangevin;

% Define all input parameters
Parameters = DefineparameterPackage(Controls);
UnpackParameters(Parameters.package,Parameters)

% Swept parameters
package = Parameters.package;
samples = unique(package(:,1));
eq_time_factors = unique(package(:,2));
N_Kuhns = unique(package(:,3));
phis = unique(package(:,4));
kds = unique(package(:,5));
Weissenbergs = unique(package(:,6));

% Define file size threshold
if BeadSpringOrMeso==0
    filesize_threshold = 1; % Modify as needed
elseif BeadSpringOrMeso==1
    filesize_threshold = 3e6; % Modify as needed
end

rerun_package = [];
% Sweep through parameters
for N_Kuhn = N_Kuhns'
 for Weissenberg = Weissenbergs'
  for phi = phis'
   for sample = samples'
    for kd = kds'
     for eq_time_factor = eq_time_factors'
         % Define filename to check
         [~,D,tau0,dtFact] = DefineTimeScale(b,length_conversion,...
             damper_conversion,BeadSpringOrMeso,phi,N_Kuhn);
         dt = tau0/dtFact;
         DefineNormalizedBindingEnergies; % (needed for filenames)
         DefineFileNames(Controls);

         if BeadSpringOrMeso==0
             filename = [output_folder,'/',OutputAtom];
         elseif BeadSpringOrMeso==1
             filename = [output_folder,'/',OutputAtom_rlx];
         end

         N_backbone_beads = Nt+(Nt-1)*(N_Kuhn-1);
         N_sidechain_beads = Nt*N_Kuhn*Ns;
         N_bead = Np*(N_backbone_beads+N_sidechain_beads);
         N_meso = Np*(Nt*Ns+Nt);

         % Check if file exists, and - if it does - whether it is large
         % enough to have successfully finished
         file_info = dir(filename);
         if ~isfile(filename) ||...
                 (isfile(filename) && file_info.bytes < filesize_threshold)
             pkg_temp = [sample eq_time_factor N_Kuhn phi kd...
                 Weissenberg N_bead N_meso];
         else
             pkg_temp = [];
         end

         rerun_package = cat(1,rerun_package,pkg_temp);
     end
    end
   end
  end
 end
end

package = rerun_package;

end
