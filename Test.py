d = {}
d.setdefault(1,[]).extend(['a','b','c'])
d.setdefault(2,[]).append('b')
count = 0
for x in d.values():
    count+=len(x)
print count
    