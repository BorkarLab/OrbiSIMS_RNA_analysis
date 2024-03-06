function parallel_assignments(ladder, spectrum, spectrum_tol, i, fragRNA, fragSize, fileID, chk, km, jm)

%constants
tol = 1E06; % tolerance 1 ppm
proton = 1.0078;
iontypes5 = [ "a-B", "a", "b", "c", "d"];
iontypes3 = ['w', 'x', 'y', 'z'];
ionloss = ["neutral", "OH", "H2O", "HPO3", "H3PO4" ];
ionlossVal = [0, 17.0027, 18.0105, 79.9663, 97.9769 ];

for l = 1:1:size(ionlossVal,2) % for type of loss - neutral, H or OH

    workSize = fragSize - ionlossVal(l);
    loss = ionloss(l);
    
    if ( chk == "segment" )
        ion = sprintf("M-%s", loss);
    end
    if ( chk == "5ion" )
        ion = sprintf("%s%d-%s", iontypes5(km),jm, loss);
    end
    if ( chk == "3ion" )
        ion = sprintf("%s%d-%s", iontypes3(km),jm, loss);
    end
    
    matchIon_parallel(workSize, spectrum, spectrum_tol, ion, loss, i, fragRNA, fragSize, fileID );
end

% ladder to maximum charge
parfor x = 1:1:ladder
    
    workSize = ( fragSize / x ) - proton;
    loss = "ladder";
    
    if ( chk == "segment" )
        ion = sprintf("M-%dH", x);
    end
    if ( chk == "5ion" )
        ion = sprintf("%s%d_%d-", iontypes5(km),jm,x);
    end
    if ( chk == "3ion" )
        ion = sprintf("%s%d_%d-", iontypes3(km),jm,x);
    end
    matchIon_parallel(workSize, spectrum, spectrum_tol, ion, loss, i, fragRNA, fragSize, fileID );
    
end
end