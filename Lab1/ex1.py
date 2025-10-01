# make an application that is able to find the alphabet of a sequence of text,this seq may be an ARN seq or ADN seq or proton seq

seq = "aabb"
entries = ""

for i in seq:
    # print(i)
    if i not in entries:
        entries = entries+str(i)

print(entries)
