# A DNA seq is given: S="ACGGGCATATGCGC" Make an app which is able to show the percentage of the components from the alphabet of the seq S.  In other words, the input of the seq S and the output is the alphabet of the seq
# and the percentage of each letter in the alphabet found in seq S

s = "ACGGGCATATGCGC"
entries = ""

for i in s:
    # print(i)
    if i not in entries:
        entries = entries+str(i)

print("\nAlphabet: " + entries)

freq = {}
for i in entries:
    freq[i] = 0

for i in s:
    freq[i] = freq[i] + 1


print("\nFrequency: ")
print(freq)

#print("Size: " + str(len(s)))

print("\nRelative freq:")
for i in range(0,len(entries)):
    percentage = round(freq[entries[i]]/len(s)*100,2)
    print(entries[i] + " " + str(percentage)+ " %")

