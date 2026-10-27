/* rx_parse.js: free-text reagent line -> stages of species + condition flags.
   "1. CH3MgBr 2. H3O+", "alc. KOH, heat", "HNO3/H2SO4", "Br2, hv", "NaNO2/HCl 0-5 C", "O3 then Zn/H2O" */
(function (root) {
'use strict';
// R groups usable inside formulas (condensed or abbreviations) -> Kekule SMILES, attach at first atom
const RG = [
  ['ch3', 'C'], ['me', 'C'], ['c2h5', 'CC'], ['ch3ch2', 'CC'], ['et', 'CC'], ['c3h7', 'CCC'], ['n-c3h7', 'CCC'], ['npr', 'CCC'], ['n-pr', 'CCC'], ['pr', 'CCC'], ['ch3ch2ch2', 'CCC'],
  ['(ch3)2ch', 'C(C)C'], ['ipr', 'C(C)C'], ['i-pr', 'C(C)C'], ['iso-c3h7', 'C(C)C'], ['i-c3h7', 'C(C)C'], ['c4h9', 'CCCC'], ['n-c4h9', 'CCCC'], ['nbu', 'CCCC'], ['n-bu', 'CCCC'], ['bu', 'CCCC'],
  ['(ch3)3c', 'C(C)(C)C'], ['tbu', 'C(C)(C)C'], ['t-bu', 'C(C)(C)C'], ['t-c4h9', 'C(C)(C)C'], ['sbu', 'C(C)CC'], ['sec-bu', 'C(C)CC'], ['ibu', 'CC(C)C'], ['i-bu', 'CC(C)C'],
  ['c6h5ch2', 'CC1=CC=CC=C1'], ['phch2', 'CC1=CC=CC=C1'], ['bn', 'CC1=CC=CC=C1'], ['benzyl', 'CC1=CC=CC=C1'], ['c6h5', 'C1=CC=CC=C1'], ['ph', 'C1=CC=CC=C1'], ['phenyl', 'C1=CC=CC=C1'],
  ['ch2=chch2', 'CC=C'], ['allyl', 'CC=C'], ['ch2=ch', 'C=C'], ['vinyl', 'C=C'], ['c6h11', 'C1CCCCC1'], ['cyclohexyl', 'C1CCCCC1'], ['r', '*'],
];
// species: id, aliases, props. Aliases are matched after compression (lowercase, no spaces/dots)
const SP = [];
const S = (id, al, p = {}) => SP.push({ id, al: al.split('|'), ...p });
// acids / water
S('HCl', 'hcl|conchcl|hydrochloricacid', { HX: 'Cl', acid: 1 }); S('HBr', 'hbr|hydrobromicacid', { HX: 'Br', acid: 1 }); S('HI', 'hi|hydroiodicacid', { HX: 'I', acid: 1 });
S('HF', 'hf'); S('H2SO4', 'h2so4|sulphuricacid|sulfuricacid', { acid: 1 }); S('H3PO4', 'h3po4|phosphoricacid', { acid: 1 }); S('HNO3', 'hno3|nitricacid|fuminghno3');
S('H3O+', 'h3o+|h3o|h+|h⁺|acid|h+/h2o|hydronium|workup', { acid: 1, workup: 1 }); S('H2O', 'h2o|water|steam|moisture'); S('D2O', 'd2o|heavywater'); S('ODm', 'od-|naod|kod|deuteroxide', { base: 1 });
S('oleum', 'oleum|fumingh2so4|so3|h2s2o7'); S('KHSO4', 'khso4', { acid: 1 }); S('pTsOH', 'ptsoh|tsoh', { acid: 1 }); S('HNO2', 'hno2|nitrousacid|nano2|nano2/hcl|kno2/hcl');
S('NaNO2', '#'); S('HOCl', 'hocl|cl2/h2o|cl2+h2o|cl2water|chlorinewater', { HOX: 'Cl' }); S('HOBr', 'hobr|br2/h2o|br2+h2o|br2water|brominewater|br2(aq)|br2aq', { HOX: 'Br' });
S('succanh', 'succinicanhydride|(ch2co)2o|butanedioicanhydride', { smi: 'O=C1CCC(=O)O1' }); S('NOCl', 'nocl|nitrosylchloride|tildenreagent|tilden'); S('IN3', 'in3|iodineazide|in₃'); S('SF4', 'sf4|dast|sulphurtetrafluoride|sulfurtetrafluoride'); S('PPh3', 'pph3|ph3p|triphenylphosphine|(c6h5)3p'); S('ClCN', 'clcn|cyanogenchloride|(cn)2|cyanogen|nccn|brcn|cyanogenbromide'); S('RS-', 'c2h5s-|ets-|ch3s-|mes-|rs-|c2h5sna|etsna|ch3sna|mesna|nasc2h5|nasch3|thiolate', { nu: 'SR' }); S('Zinin', 's2-|s--|hs-|(nh4)2s|nh4hs|na2s|na2s2|(nh4)2sx|ammoniumsulphide|ammoniumsulfide|sodiumsulphide|sodiumsulfide|zinin|h2s/nh3|h2s,nh3'); S('AlOiPr', 'al(o-tbu)3|al(otbu)3|al(o-t-bu)3|aluminiumtert-butoxide|aluminumtert-butoxide|al(o-ipr)3|al(oipr)3|al(o-i-pr)3|al(oipr₃)|aluminiumisopropoxide|aluminumisopropoxide|[(ch3)2cho]3al|al[och(ch3)2]3|mpv'); S('NO2+', 'no2+|nitronium|no2bf4|nitroniumion'); S('TsCl', 'tscl|p-tscl|tosylchloride|p-toluenesulphonylchloride|ch3c6h4so2cl'); S('SbCl5', 'sbcl5'); S('bleach', 'bleachingpowder|caocl2|ca(ocl)2|caocl₂', { hypo: 1 }); S('ArN2', 'benzenediazoniumchloride|benzenediazonium|c6h5n2cl|c6h5n2+cl-|c6h5n2+|phn2+cl-|phn2cl|phn2+|arn2+|arn2cl|diazoniumsalt'); S('chloral', 'chloral|ccl3cho|cl3ccho|trichloroethanal|trichloroacetaldehyde'); S('CuNH3', 'ammoniacalcu2cl2|ammoniacalcucl|ammoniacalcuprouschloride|cu2cl2/nh3|cu2cl2+nh3|cucl/nh3|cu(nh3)2+|[cu(nh3)2]+'); S('NH2Cl', 'nh2cl|chloramine|nh2oso3h|h2noso3h|hydroxylamine-o-sulfonicacid'); S('HCN', 'hcn|kcn/h+|nacn/h+|kcn/hcl|kcn+hcl', { cn: 1 }); S('H2O2', 'h2o2|hydrogenperoxide');
// bases / nucleophiles
S('KOH', 'koh|naoh|lioh|oh-|oh−|hydroxide|ca(oh)2|ba(oh)2|nahco3oh', { base: 1, OH: 1 });
S('NaHCO3', 'nahco3|na2co3|k2co3|bicarbonate|hco3-|hco3|co32-|co3--'); S('NaOEt', 'naoet|etona|c2h5ona|naoc2h5|sodiumethoxide|eto-|c2h5o-', { base: 1, alkoxide: 'CC' });
S('NaOMe', 'naome|meona|ch3ona|naoch3|sodiummethoxide|meo-|ch3o-', { base: 1, alkoxide: 'C' }); S('KOtBu', 'kotbu|tbuok|t-buok|(ch3)3cok|potassiumtert-butoxide|tbuo-|t-buo-', { base: 1, bulky: 1, alkoxide: 'C(C)(C)C' });
S('NaNH2', 'nanh2|sodamide|kNH2|knh2|lda|nh2-|amide-', { base: 1, strong: 1 }); S('NaH', 'nah', { base: 1 });
S('NaCN', 'nacn|kcn|cn-|cn−|cyanide', { nu: 'cn' }); S('AgCN', 'agcn', { nu: 'nc' }); S('KNO2', 'kno2', { nu: 'ono' }); S('AgNO2', 'agno2|ag+,no2-|ag+no2-|ag+/no2-|silvernitrite', { nu: 'no2' });
S('NaI', 'nai|ki|i-', { nu: 'I' }); S('AgF', 'agf|hg2f2|cof2|sbf3|naf|kf|f-|fluoride', { nu: 'F' }); S('NaSH', 'nash|khs|sh-', { nu: 'SH' }); S('NaN3', 'nan3|n3-', { nu: 'N3' });
S('NH3', 'nh3|ammonia|nh4oh|liqnh3|nh3(l)|nh3(excess)', { amine: 'N' }); S('AgNO3', 'agno3|alcagno3'); S('Ag2O', 'ag2o|moistag2o|ag2o/h2o', { OH: 1 });
S('NaOAc', 'naoac|ch3coona|ch3cook|ch3coo-|ch3cooag|agoac|acona|acok|aco-|sodiumacetate|naoocch3');
// metals / reductants
S('Na', 'na|sodium|k|potassium'); S('Mg', 'mg|magnesium'); S('Li', 'li|lithium'); S('Zn', 'zn|zinc|zndust'); S('Sn', 'sn|tin'); S('Fe', 'fe|iron'); S('Cu', 'cu|copper|cupowder');
S('MgHg', 'mg-hg|mg–hg|mg(hg)|mghg|mg/hg|magnesiumamalgam|amalgamatedmagnesium'); S('ZnHg', 'zn-hg|zn(hg)|znhg|zinc-amalgam|zincamalgam|amalgamatedzinc'); S('redP', 'redp|p(red)|p/hi|p4');
S('H2', 'h2|hydrogen|3h2|2h2'); S('Pd', 'pd|pd/c|pd-c|pt|pto2|ni|raneyni|raney-ni|rh|pt/c|catalyst', { cat: 1 }); S('Lindlar', 'lindlar|lindlarcatalyst|pd/baso4|pd-baso4|pd/caco3|pd-caco3|quinoline|pd/baso4quinoline');
S('NaBH4', 'nabh4|sodiumborohydride|lbh4|libh4', { hyd: 1 }); S('LiAlH4', 'lialh4|lah|lithiumaluminiumhydride|lithiumaluminumhydride|alh4-', { hyd: 2 }); S('DIBAL', 'dibal|dibal-h|dibalh|(ibu)2alh|i-bu2alh');
S('BH3', 'bh3|b2h6|bh3.thf|bh3thf|bh3/thf|diborane|9-bbn|sia2bh|disiamylborane'); S('NH2NH2', 'nh2nh2|n2h4|hydrazine|nh2-nh2', { amine: 'NN' }); S('SnCl2', 'sncl2|stannouschloride');
S('Hg2+', 'hgso4|hg2+|hg++|hg(oac)2|hgoac2|mercuricacetate|hg(ii)'); S('NH4Cl', 'nh4cl'); S('Na/NH3', 'na/nh3|na/liqnh3|na/nh3(l)|li/nh3|li/liqnh3|na,liqnh3|birch');
// oxidants
S('KMnO4', 'kmno4|mno4-|permanganate|baeyer|baeyersreagent'); S('K2Cr2O7', 'k2cr2o7|na2cr2o7|cr2o72-|dichromate|cro3|jones|jonesreagent|h2cro4|chromicacid'); S('PCC', 'pcc|pdc|collins|mno2');
S('OsO4', 'oso4'); S('O3', 'o3|ozone'); S('mCPBA', 'mcpba|peracid|rco3h|ch3co3h|perbenzoicacid|c6h5co3h|peroxyacid|peraceticacid|cf3co3h'); S('HIO4', 'hio4|naio4|pb(oac)4');
S('S8', 'sulphur|sulfur|s8|rhombicsulphur'); S('HCOOEt', 'hcooet|hcooc2h5|ethylformate|ethylmethanoate|hcooch2ch3', { smi: 'CCOC=O' }); S('rearr', 'rearrange|rearrangement|shift|1,2-shift|hydrideshift|methylshift|alkylshift|ringexpansion|ringcontraction|wagner-meerwein|wagnermeerwein'); S('taut', 'tautomerise|tautomerize|tautomerism|tautomerisation|tautomerization|tautomer|enolise|enolize|enolisation|enolization'); S('DMA', 'n,n-dimethylaniline|dimethylaniline|c6h5n(ch3)2|phnme2|c6h5nme2', { smi: 'CN(C)c1ccccc1' }); S('F2', 'f2|fluorine'); S('Brrad', 'br•|br·|bromineradical|brradical|brominatom'); S('couple', 'couple|coupling|dimerise|dimerize|combine|recombine|radicalcoupling|termination'); S('Schiff', 'schiff|schiffs|schiffreagent|rosaniline|fuchsin|fuchsine|magenta|decolourisedrosaniline|decolorizedfuchsin'); S('Tollens', 'tollens|tollensreagent|[ag(nh3)2]+|ag(nh3)2+|ammoniacalagno3'); S('Fehling', 'fehling|fehlings|fehlingsolution|benedict|benedicts|cu2+'); S('SeO2', 'seo2'); S('CrO2Cl2', 'cro2cl2|etard'); S('Me2S', 'me2s|(ch3)2s|dms|zn/h2o|zn,h2o|zn+h2o|zn/ch3cooh|zn/acoh|ph3p', { redwork: 1 });
S('O2', 'o2|air|oxygen'); S('aromat', 'cr2o3|v2o5|mo2o3|cr2o3/al2o3|v2o5/al2o3|mo2o3/al2o3|pt/al2o3'); S('ethene', 'ch2=ch2|ethene|ethylene|c2h4|h2c=ch2', { smi: 'C=C', dph: 1 }); S('ethyne', 'hc#ch|hc≡ch|ethyne|acetylene|c2h2', { smi: 'C#C', dph: 1 }); S('acrolein', 'ch2=chcho|acrolein|propenal|acrylaldehyde', { smi: 'C=CC=O', dph: 1 }); S('maleicanhydride', 'maleicanhydride', { smi: 'O=C1OC(=O)C=C1', dph: 1 }); S('acrylate', 'ch2=chcooch3|methylacrylate|ch2=chcooc2h5|ethylacrylate|ch2=chcooet', { smi: 'C=CC(=O)OC', dph: 1 }); S('acrylonitrile', 'ch2=chcn|acrylonitrile|vinylcyanide', { smi: 'C=CC#N', dph: 1 }); S('quinone', 'p-benzoquinone|benzoquinone|quinone', { smi: 'O=C1C=CC(=O)C=C1', dph: 1 }); S('FeTube', 'redhotirontube|redhotirontubes|redhotiron|irontube|redhotfe|redhotfetube|redhotcutube|redhottube|hotirontube');
// halogenating / others
S('Cl2', 'cl2|chlorine'); S('Br2', 'br2|bromine'); S('I2', 'i2|iodine'); S('NBS', 'nbs|n-bromosuccinimide'); S('SO2Cl2', 'so2cl2'); S('NaOCl', 'naocl|naobr|naoi|bleach|hypohalite|ocl-|clo-|obr-|bro-|oi-|hypochlorite', { hypo: 1 }); S('Cl-', 'cl-|chlorideion');
S('FeX3', 'febr3|fecl3|fe/br2|alcl3|albr3|lewisacid|anhydalcl3|anhyd.alcl3|fe', { lewis: 1 }); S('ZnCl2', 'zncl2|lucas|lucasreagent'); S('SOCl2', 'socl2|thionylchloride'); S('PCl5', 'pcl5'); S('PCl3', 'pcl3'); S('PBr3', 'pbr3|p/br2'); S('PI3', 'pi3|p/i2'); S('POCl3', 'pocl3'); S('P2O5', 'p2o5|p4o10');
S('CHCl3', 'chcl3|chloroform'); S('CCl4', 'ccl4'); S('CO2', 'co2|dryice|carbondioxide'); S('CO', 'co|co/hcl|co,hcl|gattermann-koch'); S('HCHO', 'hcho|formaldehyde|methanal|ch2o', { smi: 'C=O' });
S('CH3CHO', 'ch3cho|acetaldehyde|ethanal', { smi: 'CC=O' }); S('acetone', 'ch3coch3|(ch3)2co|acetone|propanone', { smi: 'CC(C)=O' }); S('PhCHO', 'c6h5cho|phcho|benzaldehyde', { smi: 'O=CC1=CC=CC=C1' });
S('benzene', 'c6h6|benzene', { smi: 'C1=CC=CC=C1' }); S('phenol', 'c6h5oh|phoh|phenol', { smi: 'OC1=CC=CC=C1' }); S('aniline', 'c6h5nh2|phnh2|aniline', { smi: 'NC1=CC=CC=C1', amine: 'NC1=CC=CC=C1' });
S('NH2OH', 'nh2oh|hydroxylamine|nh2oh.hcl', { amine: 'NO' }); S('PhNHNH2', 'c6h5nhnh2|phnhnh2|phenylhydrazine', { amine: 'NNC1=CC=CC=C1' }); S('DNP', '2,4-dnp|2,4dnp|dnp|2,4-dnph|bradys', { amine: 'NNC1=CC=C(C=C1[N+]([O-])=O)[N+]([O-])=O' });
S('semicarbazide', 'nh2nhconh2|nh2conhnh2|semicarbazide', { amine: 'NNC(N)=O' }); S('glycol', 'hoch2ch2oh|ethyleneglycol|ethane-1,2-diol|glycol|(ch2oh)2', { smi: 'OCCO', diol: 1 });
S('Ph3P=CH2', 'ph3p=ch2|ph3pch2|wittig|ylide', { ylide: 'C' }); S('Ph3P=CHCH3', 'ph3p=chch3', { ylide: 'CC' }); S('ethyleneoxide', 'ethyleneoxide|oxirane|epoxyethane', { smi: 'C1CO1' });
S('HBF4', 'hbf4|fluoroboricacid|nabf4'); S('H3PO2', 'h3po2|hypophosphorousacid'); S('CuCl', 'cucl|cu2cl2|cucl/hcl', { cu: 'Cl' }); S('CuBr', 'cubr|cu2br2|cubr/hbr', { cu: 'Br' }); S('CuCN', 'cucn|cucn/kcn', { cu: 'CN' });
S('naphthol', 'β-naphthol|beta-naphthol|2-naphthol|b-naphthol', { smi: 'OC1=CC=C2C=CC=CC2=C1' }); S('Hinsberg', 'c6h5so2cl|phso2cl|hinsberg|benzenesulphonylchloride|benzenesulfonylchloride');
S('HN3', 'hn3|nan3/h2so4|hydrazoicacid'); S('sodalime', 'sodalime|naoh/cao|naoh+cao|cao'); S('electrolysis', 'electrolysis|kolbeelectrolysis|electrolysed|electrolyse'); S('Al2O3', 'al2o3|alumina|tho2');
S('ROOR', 'roor|peroxide|(phcoo)2|(c6h5coo)2o2|benzoylperoxide|bz2o2|(ch3)3cooc(ch3)3|h2o2/hbr|kharasch');
S('NaHSO3', 'nahso3|sodiumbisulphite|sodiumbisulfite'); S('Br2/NaOH', 'br2/naoh|br2/koh|br2+naoh|br2+koh|naobr|kobr|hofmann', { hofmann: 1 });
S('Zn/HCl', 'zn/hcl|zn+hcl|zn,hcl'); S('Sn/HCl', 'sn/hcl|sn+hcl|fe/hcl|fe+hcl|sn,hcl|fe,hcl|zn/hcl'); S('Zn/NH4Cl', 'zn/nh4cl'); S('Zn/NaOH', 'zn/naoh');
S('alcKOH', 'alc.koh|alckoh|alcoholickoh|ethanolickoh|koh/etoh|koh/c2h5oh|koh(alc)|koh(alc.)|alc.naoh|alcnaoh|koh,ethanol|naoh/etoh', { base: 1, OH: 1, alc: 1 });
S('aqKOH', 'aq.koh|aqkoh|aqueouskoh|aq.naoh|aqnaoh|koh(aq)|naoh(aq)|aqueousnaoh|koh/h2o|naoh/h2o|dil.naoh|dilnaoh|dil.koh|dilkoh', { base: 1, OH: 1, aq: 1 });
S('Ac2O', '(ch3co)2o|ac2o|aceticanhydride|(ch3co)2o/ch3coona', { acyl: 'CC(=O)', smi: 'CC(=O)OC(C)=O' });
S('AcCl', 'ch3cocl|accl|acetylchloride|ch3cocl/pyridine|ch3coi|acetyliodide|ch3cobr|acetylbromide', { acyl: 'CC(=O)', smi: 'CC(=O)Cl' }); S('BzCl', 'c6h5cocl|phcocl|benzoylchloride|bzcl', { acyl: 'C(=O)C1=CC=CC=C1', smi: 'O=C(Cl)C1=CC=CC=C1' });
S('AcOH', 'ch3cooh|aceticacid|acoh|glacialaceticacid|hcooh', { smi: 'CC(O)=O' });
// flags
const FLAGS = { heat: 'heat|δ|∆|delta|reflux|warm|boil|hot|hightemp|hightemperature|hight|fuse|fusion|pyrolysis', light: 'hv|hν|hu|light|uv|sunlight|diffusedsunlight|irradiation|lightenergy',
  cold: 'cold|ice|ice-cold|lowtemp|lowtemperature|lowt|0°c|0c|-78c|-78', dil: 'dil|dilute', conc: 'conc|concentrated|concd', aq: 'aq|aqueous', alc: 'alc|alcoholic|ethanolic|etoh|c2h5oh|ethanol|meoh|ch3oh|methanol|isopropanol|propan-2-ol|2-propanol|iproh|ipa',
  dry: 'dry|anhydrous|anhyd', ether: 'ether|et2o|dryether|diethylether|thf|dioxane', excess: 'excess|xs|exc', one: '1eq|1mole|1mol|oneequivalent|1equiv', pyr: 'pyridine|py|et3n', acetone: 'acetone',
  ccl4: 'ccl4|ch2cl2|cs2|chcl3solvent|inert|xylene|dryxylene', polar: 'dmso|dmf|hmpa', pressure: 'pressure|highpressure|300atm|atm|bar', glycol: 'ethyleneglycol|glycol|hoch2ch2oh|(ch2oh)2', red: 'red', liq: 'liq.|liq|liquid|liquefied', dark: 'dark|indark|absenceoflight|nolight', cat: '(cat.)|(cat)|cat.|catalyst|catalytic|trace|traces', basecat: 'base(cat.)|base(cat)|basecat.|basecatalysed|basecatalyzed', protic: 'polarprotic|proticsolvent', aprotic: 'polaraprotic|aproticsolvent', moist: 'moist|moisture|wet', sealed: 'sealedtube|sealed', rt: 'rt|r.t.|roomtemperature|roomtemp', lewis: 'lewisacid' };
const ALIAS = new Map();
for (const s of SP) for (const a of s.al) if (a && a !== '#') ALIAS.set(a, s.smi ? { sp: s.id, smi: s.smi } : { sp: s.id });
for (const k in FLAGS) for (const a of FLAGS[k].split('|')) if (!ALIAS.has(a)) ALIAS.set(a, { flag: k });
// formula families with R groups (named species win)
const setA = (k, v) => { if (!ALIAS.has(k)) ALIAS.set(k, v); };
for (const [r, smi] of RG) {
  for (const x of ['cl', 'br', 'i']) {
    const X = x === 'cl' ? 'Cl' : x === 'br' ? 'Br' : 'I';
    setA(r + 'mg' + x, { sp: 'RMgX', smi: smi === '*' ? 'C' : smi, X, gen: smi === '*' }); setA(r + x, { sp: 'RX', smi: (smi === '*' ? 'C' : smi), X, gen: smi === '*' });
  }
  setA(r + 'mgx', { sp: 'RMgX', smi: smi === '*' ? 'C' : smi, X: 'Br', gen: smi === '*' }); setA(r + 'x', { sp: 'RX', smi: smi === '*' ? 'C' : smi, X: 'Cl', gen: smi === '*' });
  setA(r + 'li', { sp: 'RLi', smi: smi === '*' ? 'C' : smi });
  setA(r + 'oh', { sp: 'ROH', smi: smi === '*' ? 'C' : smi }); setA(r + 'nh2', { sp: 'RNH2', smi: smi === '*' ? 'C' : smi });
  if (smi !== '*') { setA('na' + 'o' + r, { sp: 'NaOR', smi }); setA(r + 'ona', { sp: 'NaOR', smi }); setA(r + 'ok', { sp: 'NaOR', smi }); setA(r + 'cocl', { sp: 'RCOCl', smi }); setA('(' + r + ')2culi', { sp: 'R2CuLi', smi }); setA('(' + r + ')2cd', { sp: 'R2Cd', smi }); setA(r + '2cd', { sp: 'R2Cd', smi }); setA(r + '2culi', { sp: 'R2CuLi', smi }); }
}
ALIAS.set('ch3i', { sp: 'RX', smi: 'C', X: 'I' }); ALIAS.set('mei', { sp: 'RX', smi: 'C', X: 'I' }); ALIAS.set('ch3oh', { sp: 'ROH', smi: 'C' }); ALIAS.set('etoh', { sp: 'ROH', smi: 'CC' }); ALIAS.set('c2h5oh', { sp: 'ROH', smi: 'CC' });
ALIAS.set('(ch3)2nh', { sp: 'R2NH', smi: 'C' }); ALIAS.set('me2nh', { sp: 'R2NH', smi: 'C' }); ALIAS.set('(c2h5)2nh', { sp: 'R2NH', smi: 'CC' });
ALIAS.set('ch2n2', { sp: 'CH2N2' }); ALIAS.set('r2cd', { sp: 'R2Cd', smi: 'C' }); ALIAS.set('cdcl2', { sp: 'CdCl2' }); ALIAS.set('diazomethane', { sp: 'CH2N2' });
const KEYS = [...ALIAS.keys()].sort((a, b) => b.length - a.length);
const MAXK = KEYS[0].length;

function normalize(t) {
  const SUB = '₀₁₂₃₄₅₆₇₈₉', SUP = '⁰¹²³⁴⁵⁶⁷⁸⁹';
  let s = '';
  for (const ch of t) { let k = SUB.indexOf(ch); if (k >= 0) { s += k; continue; } k = SUP.indexOf(ch); if (k >= 0) { s += k; continue; } s += ch; }
  return s.replace(/⁺/g, '+').replace(/[⁻−–—]/g, '-').replace(/[→⟶➝]/g, ' then ').replace(/·/g, '.').replace(/\u00a0/g, ' ');
}
function splitStages(t) {
  t = normalize(t);
  const parts = t.split(/\bthen\b|;|\n|\bfollowed by\b|\(\s*(?:i{1,3}|iv|v)\s*\)|(?:^|\s)(?:i{1,3}|iv)\)|(?:^|[\s,])\d\s*[.)](?!\d)/i).map(x => x && x.trim()).filter(x => x);
  return parts.length ? parts : [''];
}
function temps(s) {
  const out = []; const re = /(-?\d{1,4})\s*(?:-\s*(-?\d{1,4})\s*)?(°\s*c|°c|deg\s*c|k\b|°|c\b)/gi; let m;
  while ((m = re.exec(s))) { let v = +m[1], v2 = m[2] != null ? +m[2] : v; const u = m[3].toLowerCase().replace(/\s/g, ''); let K = u === 'k' ? (v + v2) / 2 : (v + v2) / 2 + 273; out.push({ K, txt: m[0] }); }
  return out;
}
function parseStage(txt) {
  const st = { raw: txt, sp: new Map(), flags: new Set(), T: null, unknown: [], eq: null };
  let s = txt.toLowerCase();
  // phrases first (they contain little words that are separators elsewhere)
  s = s.replace(/\bin\s+(?:the\s+)?dark\b|\babsence\s+of\s+light\b|\bno\s+light\b/g, ' dark ').replace(/\b(?:in\s+)?(?:the\s+)?presence\s+of\b/g, ' , ')
       .replace(/\broom\s+temp(?:erature)?\b/g, ' rt ').replace(/\bsealed\s+tube\b/g, ' sealed ').replace(/\bred[\s-]*hot\s+(iron|fe|cu|copper)\s*(tubes?)?\b/g, ' redhotirontube ');
  // amounts: "2 eq", "(1 equiv.)", "2 mol", "excess", "xs", leading coefficients "2HBr" / "3 H2"
  s = s.replace(/\(?\s*(\d+(?:\.\d+)?)\s*(?:eq|equiv|equivalents?|moles?|mol)\b\.?\s*\)?/g, (m0, n) => { st.eq = +n; return ','; });
  s = s.replace(/\(?\b(?:excess|xs|in\s+excess|large\s+excess)\b\)?/g, () => { st.eq = Infinity; return ',excess,'; });
  s = s.replace(/(^|[\s,(+\/])([2-6])\s*(?=(?:h2|hbr|hcl|hi|br2|cl2|i2|hno3|naoh|koh|nanh2|ch3i|mei|ch3mgbr|ch3mgi|rmgx|o3|nabh4|lialh4|na|naome|naoet|ch3cl|ch3br)\b)/g, (m0, pre, n) => { if (st.eq == null) st.eq = +n; return pre; });
  s = s.replace(/\d+(?:\.\d+)?\s*%/g, ' ');
  s = s.replace(/(\d+(?:\s*-\s*\d+)?)\s*(?:atm|bar)\b/g, () => { st.flags.add('pressure'); return ' '; });
  // apostrophes and filler words: "Schiff's reagent", "Tollens' reagent", "Fehling's solution", "Lucas test"
  s = s.replace(/['’‘`]/g, '').replace(/\b(?:reagents?|solution|soln\.?|test)\b/g, ' ');
  // little linking words are separators, not reagents ("Na in dry ether", "Br2 with FeBr3")
  s = s.replace(/\b(?:in|with|and|using|under|over|catalyst of|by|on|at)\b/g, ' , ');
  for (const t of temps(s)) { st.T = st.T == null ? t.K : st.T; s = s.replace(t.txt.toLowerCase(), ' '); }
  if (/\b(\d{3})\s*$/.test(s) && st.T == null) { st.T = +RegExp.$1; s = s.replace(/\b\d{3}\s*$/, ''); }
  const c = s.replace(/\s+/g, '').replace(/,(?=\d)/g, '#').replace(/\.(?![a-z]?\))/g, '.');
  // greedy tokenization; separators , / + : are boundaries but aliases may contain them
  let i = 0, junk = '';
  const flushJunk = () => { const j = junk.replace(/[\s,/+.:()\-]/g, ''); if (j) st.unknown.push(junk.replace(/^[,/+.:]+|[,/+.:]+$/g, '')); junk = ''; };
  const cc = c.replace(/#/g, ',');
  while (i < cc.length) {
    let hit = null;
    for (let L = Math.min(MAXK, cc.length - i); L > 0; L--) { const sub = cc.substr(i, L); const v = ALIAS.get(sub); if (v) { // word boundary check for very short aliases
        if (L <= 2 && /[a-z0-9]/.test(cc[i + L] || '') && !/^(h\+|oh|na|k|zn|sn|fe|cu|mg|li|pd|ni|pt)$/.test(sub)) { const nx = cc[i + L]; if (/[a-z]/.test(nx)) continue; }
        hit = [sub, v]; break; } }
    if (!hit) { junk += cc[i]; i++; continue; }
    flushJunk();
    const [sub, v] = hit; i += sub.length;
    if (v.flag) st.flags.add(v.flag);
    else {
      const prev = st.sp.get(v.sp); const rec = { ...v, txt: sub };
      if (!prev) st.sp.set(v.sp, rec); else if (!Array.isArray(prev)) st.sp.set(v.sp, [prev, rec]); else prev.push(rec);
    }
  }
  flushJunk();
  if (st.flags.has('one') && st.eq == null) st.eq = 1; if (st.flags.has('excess') && st.eq == null) st.eq = Infinity;
  // derived combos
  const has = k => st.sp.has(k);
  if (has('H2SO4') && has('HNO3')) st.sp.set('NO2+', { sp: 'NO2+' });
  if (has('NaNO2') || has('HNO2')) st.sp.set('HNO2', { sp: 'HNO2' });
  if ((has('Zn') || has('Sn') || has('Fe')) && has('HCl')) st.sp.set('Sn/HCl', { sp: 'Sn/HCl' });
  if (has('Zn') && has('H2O')) st.sp.set('Me2S', { sp: 'Me2S' });
  if (has('Na') && (has('NH3') || st.flags.has('red'))) { if (has('NH3')) st.sp.set('Na/NH3', { sp: 'Na/NH3' }); }
  if (has('Br2') && has('KOH') && !st.flags.has('alc')) st.sp.set('Br2/NaOH', { sp: 'Br2/NaOH' });
  if (has('I2') && (has('KOH') || has('NaHCO3'))) st.sp.set('NaOCl', { sp: 'NaOCl', X: 'I' });
  if (has('Cl2') && has('KOH')) st.sp.set('NaOCl', { sp: 'NaOCl', X: 'Cl' });
  if (has('bleach')) st.sp.set('NaOCl', { sp: 'NaOCl', X: 'Cl' });
  if (has('ZnCl2') && /lucas/i.test(st.raw)) st.sp.set('HCl', { sp: 'HCl' });
  if (has('HCl') && has('ZnCl2')) st.sp.set('Lucas', { sp: 'Lucas' });
  if (/\bdil/i.test(st.raw)) st.flags.add('dil'); if (/\bconc/i.test(st.raw)) st.flags.add('conc');
  if (has('Br2') && st.flags.has('red')) st.sp.set('HVZ', { sp: 'HVZ' });
  if (has('Br2') && has('redP')) st.sp.set('HVZ', { sp: 'HVZ' });
  if (has('H2O2') && has('KOH')) st.sp.set('H2O2/OH', { sp: 'H2O2/OH' });
  if (has('KOH') && (st.flags.has('alc') || (has('ROH') && !has('H2O') && !st.flags.has('aq')))) st.sp.set('alcKOH', { sp: 'alcKOH' });
  if (has('KOH') && (st.flags.has('aq') || has('H2O'))) st.sp.set('aqKOH', { sp: 'aqKOH' });
  if (has('KOH') && st.flags.has('dil')) st.sp.set('aqKOH', { sp: 'aqKOH', dil: 1 });
  if (has('KMnO4') && (st.flags.has('cold') || st.flags.has('dil'))) st.flags.add('mild');
  if (has('NaH') || has('Na')) st.flags.add('metal');
  if (has('Pd') && has('H2') === false && has('Lindlar') === false) { /* Pd alone */ }
  return st;
}
/* a comma list that puts a protic / acidic reagent BEFORE an organometallic cannot all be in one flask: split it there */
const ORGM = ['RMgX', 'RLi', 'LiAlH4', 'R2CuLi', 'DIBAL'], PROT = ['H3O+', 'H2O', 'glycol', 'H2SO4', 'HCl', 'TsOH', 'ROH'];
function orderSplit(txt) {
  const tok = txt.split(',').map(x => x.trim()).filter(Boolean); if (tok.length < 2) return [txt];
  const sps = tok.map(t => { try { return parseStage(t).sp; } catch (e) { return new Map(); } });
  const i = sps.findIndex(sp => ORGM.some(k => sp.has(k))); if (i <= 0) return [txt];
  if (!sps.slice(0, i).some(sp => PROT.some(k => sp.has(k)))) return [txt];
  return [tok.slice(0, i).join(', '), tok.slice(i).join(', ')];
}
function parse(text) {
  const stages = splitStages(text || '').flatMap(orderSplit).map(parseStage);
  return { stages, unknown: stages.flatMap(s => s.unknown) };
}
function stageCtx(st) {
  const f = st.flags;
  return {
    st, raw: st.raw,
    has: k => st.sp.has(k), any: ks => ks.some(k => st.sp.has(k)), get: k => { const v = st.sp.get(k); return Array.isArray(v) ? v[0] : v; },
    flag: k => f.has(k), T: st.T, eq: st.eq,
    heat: f.has('heat') || (st.T != null && st.T >= 353), light: f.has('light'), cold: f.has('cold') || (st.T != null && st.T <= 283),
  };
}
root.RxParse = { parse, stageCtx, normalize, splitStages, ALIAS, RG, SP };
})(typeof window !== 'undefined' ? window : globalThis);
