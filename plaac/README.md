# PLAAC

Not committed — build it here. There are no releases and no prebuilt jar.

```bash
git clone --depth 1 https://github.com/whitehead/plaac.git /tmp/plaac
cp -r /tmp/plaac/cli/src plaac/
mkdir -p plaac/example && cp /tmp/plaac/cli/example/four_classic_prions.fasta plaac/example/
cp -r /tmp/plaac/web/bg_freqs plaac/

cd plaac && mkdir -p target
javac -source 8 -target 8 -nowarn -d target src/*.java
cp src/mainClass target/. && cp -r src/util target/.
cd target && jar cmf mainClass plaac.jar * && cd ..
mv target/plaac.jar plaac.jar && rm -rf target
```

Needs a JDK; `openjdk` is in `environment.yml`.

**Check it works before trusting any output.** The four classic prions should
score LLR 29–51 with a domain called for each:

```bash
java -jar plaac/plaac.jar -i plaac/example/four_classic_prions.fasta -a 1.0 -c 60
```

Expected: Sup35p 51.2, Rnq1p 47.5, Mot3p 40.8, Ure2p 29.2.

Cite Lancaster et al. (2014) *Bioinformatics* and Alberti et al. (2009) *Cell*.
