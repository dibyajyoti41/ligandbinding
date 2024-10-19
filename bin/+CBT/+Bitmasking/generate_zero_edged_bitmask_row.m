function bitmask = generate_zero_edged_bitmask_row(bitmaskLength, edgeLength)
    bitmask = true(1, bitmaskLength);
    bitmask([1:min(edgeLength, bitmaskLength), max(bitmaskLength + 1 - edgeLength, 1):bitmaskLength]) = false;
end