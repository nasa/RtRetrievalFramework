------------------------------------------------------------
--- Helper functions that add additional non L2 specific
--- functionality
------------------------------------------------------------

function table.contains(table, element)
  for _, value in pairs(table) do
    if value == element then
      return true
    end
  end
  return false
end

function table.index(table, element)
  for idx, value in pairs(table) do
    if value == element then
      return idx
    end
  end
  return -1
end
