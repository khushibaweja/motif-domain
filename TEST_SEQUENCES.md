# Test Sequences for Real API

## Quick Test Sequences

### 1. Short Protein (Fast - ~30 seconds)

**Human Insulin B Chain**

```
FVNQHLCGSHLVEALYLVCGERGFFYTPKA
```

Expected: Signal peptide, insulin family signature

### 2. Medium Protein (Medium - ~60 seconds)

**P53 DNA Binding Domain Fragment**

```
SQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD
```

Expected: P53 DNA binding domain, P53 tetramerization domain

### 3. Known Kinase (Medium - ~60 seconds)

**Protein Kinase Fragment**

```
MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREISLLKELKDDNIVRLYDIVHSDAHKLYLVFEFLDLDLKRYMEGIPKDQPLGADIVKKFMMQLCKGIAYCHSHRILHRDLKPQNLLINTEGAIKLADFGLAIEVQGDQPSSTVRRVSLPKQEKMSAKSGPLYLGTEYIYNEYVEQVGQHQRSQG
```

Expected: Protein kinase domain, ATP binding site

### 4. Transcription Factor (Medium - ~60 seconds)

**Zinc Finger Protein Fragment**

```
MGKVRFKCYTCGKFGQSSHLKAHMKHEGRTHTGEKPYRCNECGKAFIQSSNLKQHQKTHTGEKPFQCRICDKTFSRSGALQAHIKRRHGKTSKCPDCAKNFHQTKNLIQHRKIHTKEKP
```

Expected: Multiple C2H2 zinc fingers

### 5. From NCBI - Try These Accessions

**Human Hemoglobin Alpha (Fast)**

```
Accession: P69905
Expected: Globin domain
Time: ~45 seconds
```

**Human Growth Hormone (Fast)**

```
Accession: P01241
Expected: Growth hormone family signature
Time: ~40 seconds
```

**Epidermal Growth Factor Receptor (Slow)**

```
Accession: P00533
Expected: Multiple domains (kinase, EGF-like, transmembrane)
Time: ~90 seconds
```

**Green Fluorescent Protein (Fast)**

```
Accession: P42212
Expected: GFP-like domain
Time: ~35 seconds
```

## Testing Checklist

- [ ] Test with short sequence (~30 aa)
- [ ] Test with medium sequence (~200 aa)
- [ ] Test with NCBI accession number
- [ ] Test with DNA sequence (will be translated)
- [ ] Test with FASTA file upload
- [ ] Verify domains appear in visualization
- [ ] Check hover tooltips work on interactive plot
- [ ] Verify domain information shows database source
- [ ] Test with sequence that has no domains
- [ ] Test timeout with very long sequence

## Expected Results

### Good Results:

✅ Multiple domains detected
✅ Domains shown in both visualizations  
✅ Hover shows domain details
✅ Database source indicated (Pfam, SMART, etc.)
✅ Accurate positions
✅ Success message with domain count

### No Domains Found (Normal):

✅ Message: "No domains detected"
✅ Protein features still shown
✅ No visualizations displayed
✅ This is expected for some sequences!

### Error Cases:

❌ Timeout after 3 minutes
❌ Network error (check internet/server)
❌ Invalid sequence characters
❌ Sequence too long (>40,000 aa)

## Notes

- **First analysis**: May take longer as EBI servers initialize
- **Subsequent analyses**: Usually faster due to caching
- **Console output**: Check terminal for detailed progress
- **Rate limiting**: EBI may limit if many requests in short time

## Pro Tips

1. **Start small**: Test with short sequences first
2. **Check console**: Watch Flask terminal for progress
3. **Be patient**: Real analysis takes 30-120 seconds
4. **Known proteins**: Use NCBI accessions for proteins with known domains
5. **Browser console**: Open DevTools to see any JavaScript errors

Happy testing! 🧬
