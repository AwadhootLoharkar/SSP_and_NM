function [reflection, transmission, absorption] = RTA_iter(energyStep, energy, damping, potentialPositions, potentialValues, delta, energyScale)
    % Calculate reflection, transmission, and absorption coefficients
    [~, numPos] = size(potentialPositions); 
    waveNumInit = sqrt((energy + 1i * damping) / energyScale); 
    waveNumFinal = sqrt((energy + 1i * damping - energyStep) / energyScale); 
    transCoeff = 2 * waveNumInit / (waveNumInit + waveNumFinal); 
    reflCoeff = (waveNumInit - waveNumFinal) / (waveNumFinal + waveNumInit); 
    potMatrix = zeros(numPos, numPos); 
    
    for ii = 1:numPos
        waveFuncRight(ii) = transCoeff * exp(1i * waveNumFinal * potentialPositions(ii)); 
        potMatrix(ii, ii) = delta * potentialValues(ii) / energyScale; 
        
        for jj = 1:numPos
            greensFunc(ii, jj) = GreensFun(energyStep, potentialPositions(ii), potentialPositions(jj), energy, damping, energyScale); 
        end
    end
    
    scatterMatrix = eye(numPos) - greensFunc * potMatrix; 
    waveFuncSol = scatterMatrix \ (waveFuncRight.'); 
    reflWaveFunc = reflCoeff * exp(-1i * waveNumInit * (potentialPositions(1) - 1)) + ExtraTerm(potentialPositions(1) - 1, potentialPositions, potentialValues, delta, waveFuncSol, energyStep, energy, damping, energyScale); 
    transWaveFunc = transCoeff * exp(1i * waveNumFinal * (potentialPositions(numPos) + 1)) + ExtraTerm(potentialPositions(numPos) + 1, potentialPositions, potentialValues, delta, waveFuncSol, energyStep, energy, damping, energyScale); 
    
    transmission = (real(waveNumFinal) / real(waveNumInit)) * abs(transWaveFunc)^2; 
    reflection = (abs(reflWaveFunc / (exp(1i * waveNumInit * (potentialPositions(1) - 1)))))^2; 
    absorption = 1 - (reflection + transmission); 
end