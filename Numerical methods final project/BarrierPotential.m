function barrierPotential = BarrierPotential(positions, barrierStart, barrierEnd, energy)
    % Compute the barrier potential for given positions
    barrierPotential = zeros(size(positions));
    barrierPotential(barrierStart <= positions & positions <= barrierEnd) = energy;
end