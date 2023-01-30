refs = snakemake.params["refs"]
primers = snakemake.params["primers"]

if os.path.exists(refs):
    targets = os.listdir(readsDir)
    for target in targets:
        print(target)
else:
    sys.exit("\t**Reads not found**")