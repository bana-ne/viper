f = open(options.gtf)
i = 0;
for l in f:
    if l.startswith("#"):
        continue
    else:
        if i < 1:
            i+=1
            tmp = l.strip().split("\t")
            #select out last col -- gene attributes list
            attribs = tmp[-1].strip().split(";")
            #find the gene id attrib
            refID = list(filter(lambda x: x.startswith('gene_id'), attribs))
            print(refID)
        if refID:
            refID = refID[0].split(" ")[1]
            #NOTE: we get something like "NM_001195025"; need to eval as
            refID = eval(refID)
            
            #DROP the old gene_id attrib
            attribs = [a for a in attribs if not a.startswith("gene_id")]
            #INSERT the new gene_id at the head of the list
            attribs.insert(0, "gene_id \"%s\"" % (dict[refID]))
            
            #INSERT back into tmp
            tmp[-1] = ";".join(attribs)

        print("\t".join(tmp))
        
f.close()

f = open(options.gtf)
for l in f:
    if l.startswith("#"):
        continue
    else:
        tmp = l.strip().split("\t")
        #select out last col -- gene attributes list
        attribs = tmp[-1].strip().split(";")
        #find the gene id attrib
        refID = filter(lambda x: x.startswith('gene_id'), attribs)
        if refID:
            refID = refID[0].split(" ")[1]
            #NOTE: we get something like "NM_001195025"; need to eval as
            refID = eval(refID)
            
            #DROP the old gene_id attrib
            attribs = [a for a in attribs if not a.startswith("gene_id")]
            #INSERT the new gene_id at the head of the list
            attribs.insert(0, "gene_id \"%s\"" % (dict[refID]))
            
            #INSERT back into tmp
            tmp[-1] = ";".join(attribs)

        print("\t".join(tmp))

f.close()

