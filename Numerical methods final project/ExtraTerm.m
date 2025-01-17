function [extraTermSolution] = ExtraTerm(position, potentialPositions, potentialValues, delta, waveFunctionSolution, energyStep, energy, damping, energyScale)
    % Calculate the extra term in the wave function solution
    extraTermSolution = 0;
    [~, numPotentialPositions] = size(potentialPositions);
    scale = delta / energyScale;

    % Iterate over each potential position
    parfor i = 1:numPotentialPositions
        extraTermSolution = extraTermSolution + GreensFun(energyStep, position, potentialPositions(i), energy, damping, energyScale) * scale * potentialValues(i) * waveFunctionSolution(i);
    end
end