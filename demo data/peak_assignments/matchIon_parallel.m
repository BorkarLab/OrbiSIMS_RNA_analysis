function matchIon_parallel(workSize, spectrum, spectrum_tol, ion, loss, i, fragRNA, fragSize, fileID )

tol = 1E06; % 1 ppm
temp = nan(size(spectrum,1));

% tolwork = workSize / tol * 2; % 2 ppm error
temp = (abs(spectrum - workSize ) <= spectrum_tol);


if ( sum(temp) > 0 )
    matrix = spectrum(temp>0);
     for k = 1:1:sum(temp(:)) % if more than one match
        str = [ num2str(i), fragRNA, num2str(fragSize), num2str(workSize), matrix(k,1) ion loss ];
        writematrix(str,fileID,'WriteMode','append', 'Delimiter','tab');
    end
end

