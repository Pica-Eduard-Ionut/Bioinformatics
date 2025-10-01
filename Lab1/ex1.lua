local seq = "aabb"
local entries = ""

for i = 1, #seq do
    local char = seq:sub(i, i)
    if not entries:find(char, 1, true) then
        entries = entries .. char
    end
end

print(entries)
