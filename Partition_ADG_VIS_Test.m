%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Test script for the ADG extrapolation model code. The ADG ext code is run
%for one specified input of ag and adg at specified hyperspectral
%wavelengths (lambda). The resulting output from the test script is saved
%to ADG_part_test_run_yyyymmdd.xls for comparison with the provided output
%file ADG_part_test_run.xls.
%
%Reference:
%
%Kehrli, M. D., Stramski, D., Reynolds, R. A., & Joshi, I. D. (2023).
%Estimation of chromophoric dissolved organic matter and non-algal
%particulate absorption coefficients of seawater in the ultraviolet by
%extrapolation from the visible spectral region. Optics Express, 31(11),
%17450. https://doi.org/10.1364/OE.486354
%
%Created: July 31, 2023
%Completed: August 14, 2023
%Updates: 
%
%M. D. Kehrli, D. Stramski, R. A. Reynolds, and I. D. Joshi
%Ocean Optics Research Laboratory, Scripps Institution of Oceanography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear command window and workspace; close figures
clc; clearvars; close all;

%define input parameters

%input adg [m^-1]
adg = [0.191962469768812	0.189172411984555...
    0.186366843375226	0.183531206057442	0.180686560246946...
    0.177775056783389	0.175624444838783	0.173124652582223...
    0.170661172369395	0.168234017424789	0.165843428477267...
    0.163489707555984	0.161172379445197	0.158890370421873...
    0.156642219010794	0.154426349194468	0.152241479669686...
    0.150086158853360	0.147958894276100	0.145858681748404...
    0.143785518634701	0.141739502719249	0.139721035777558...
    0.137731033561815	0.135770263325663	0.133839092701756...
    0.131937697695926	0.130065149716334	0.128219864112835...
    0.126400433372310	0.124604960803797	0.122832053748404...
    0.121081141914584	0.119351898684264	0.117644734019541...
    0.115960794130328	0.114300928853360	0.112666365754235...
    0.111058122281931	0.109476891459774	0.107923470561815...
    0.106398375241115	0.104901577287762	0.103433246535576...
    0.101992850768812	0.100579622206130	0.0991928849466544...
    0.0978314598008818	0.0964940297746428	0.0951796976696865...
    0.0938871690486952	0.0926154305880539	0.0913639318679373...
    0.0901317441274124	0.0889185219058381	0.0877242791682288...
    0.0865487080982579	0.0853919610749343	0.0842545033810568...
    0.0831361127629810	0.0820367216871792	0.0809560519583162...
    0.0798932948125437	0.0788478664568585	0.0778188881390743...
    0.0768051511828060	0.0758059467688119	0.0748208928591909...
    0.0738492959320772	0.0728911909670626	0.0719466753344095...
    0.0710158854422812	0.0700992833985495	0.0691973337892200...
    0.0683099252965086	0.0674371584743512	0.0665785888650218...
    0.0657333240836807	0.0649006203110859	0.0640793965355757...
    0.0632685636667710	0.0624674567513191	0.0616755282673541...
    0.0608922028416982	0.0601173802265378	0.0593507271886369...
    0.0585920914014649	0.0578411653519022	0.0570977109379081...
    0.0563616035530684	0.0556332016288702	0.0549128979087536...
    0.0542014339058381	0.0534996908387827	0.0528084141798906...
    0.0521281302556923	0.0514589601449052	0.0508006665559839...
    0.0501527633431559	0.0495146287105028	0.0488855063664795...
    0.0482651699145845	0.0476533869816398	0.0470501266755174...
    0.0464554634131267	0.0458693714160422	0.0452915778212900...
    0.0447217698037973	0.0441593028562754	0.0436037780137098...
    0.0430551116084620	0.0425131681828060	0.0419782052615232...
    0.0414507807454882	0.0409313507863046	0.0404201543664795...
    0.0399172659554008	0.0394222846259547	0.0389348703023396...
    0.0384544427279955	0.0379803981361588	0.0375122635909693...
    0.0370499045997157	0.0365931287163337	0.0361419560195407...
    0.0356963054451967	0.0352560857658964	0.0348210783081705...
    0.0343912011244970	0.0339663898212900	0.0335466637221646...
    0.0331321479087535	0.0327231872877623	0.0323201003198323...
    0.0319232003431559	0.0315327168562754	0.0311487488912608...
    0.0307713970749343	0.0304005525472375	0.0300360171332433...
    0.0296778726142929	0.0293263582994241	0.0289816433519022...
    0.0286438625938848	0.0283130698621063	0.0279893056026311...
    0.0276723057542346	0.0273614792148760	0.0270562570399489...
    0.0267560041361588	0.0264601288300364	0.0261681915414066...
    0.0258798410078789	0.0255950950632725	0.0253140259466544...
    0.0250363097600655	0.0247619418271209	0.0244909469699780...
    0.0242232086434474	0.0239587624364503	0.0236976595676457...
    0.0234398934043804	0.0231857908533600	0.0229352968970918...
    0.0226883872556923	0.0224451799233308	0.0222057192586078...
    0.0219697429145845	0.0217373611507361	0.0215086816959256...
    0.0212838857484037	0.0210627661798906	0.0208451329583162...
    0.0206304734160422	0.0204186380370334	0.0202092088533600...
    0.0200019860166253	0.0197969776113775	0.0195948874335349...
    0.0193962164772667	0.0192018511303279	0.0190123726347011...
    0.0188284560078789	0.0186502122411151	0.0184774895297448...
    0.0183096100282871	0.0181460536784329	0.0179861526317856...
    0.0178294225968002	0.0176752003373250	0.0175232132148760...
    0.0173736675997157	0.0172272145209985	0.0170842339845553...
    0.0169448525414066	0.0168089312994241	0.0166763140107944...
    0.0165464782352842	0.0164183302411151	0.0162902592819314...
    0.0161610231244970	0.0160300177192492	0.0158971226259547...
    0.0157625736113775	0.0156272391157506	0.0154922949524853...
    0.0153590662877623	0.0152287000516107	0.0151018083752259...
    0.0149784647804737	0.0148583278387827	0.0147405559641471...
    0.0146240971070043	0.0145082282265378	0.0143923761886369...
    0.0142763521361588	0.0141603407134183	0.0140448231070043...
    0.0139302751419897	0.0138174041099197	0.0137064998504445...
    0.0135977173956340	0.0134909102469460	0.0133858559379081...
    0.0132823141857215	0.0131804351186661	0.0130803496026311...
    0.0129824620137098	0.0128871599029226	0.0127948790399489...
    0.0127058510924270	0.0126200168387827	0.0125369586551092...
    0.0124561313460713	0.0123767490341180	0.0122980414510276...
    0.0122194256988410	0.0121405942527769	0.0120616606026311...
    0.0119827954830976	0.0119042493577331	0.0118265638387827...
    0.0117502784043804	0.0116754878737681	0.0116020600370334...
    0.0115294517192492	0.0114572780545261	0.0113850607688119...
    0.0113120050107944	0.0112373858387827	0.0111610816726020...
    0.0110831480166253	0.0110040121040888	0.0109240616755174...
    0.0108440694568585	0.0107651710428643	0.0106883580574416...
    0.0106142769466544	0.0105436244568585	0.0104769349554008...
    0.0104145802207069	0.0103564408504445	0.0103020829174999...
    0.0102513720049635	0.0102039601915524	0.0101589754102113...
    0.0101152440312025	0.0100714710282871	0.0100261869262463...
    0.00997794272216463	0.00992497960263110	0.00986618510117338...
    0.00980135814198970	0.00973099777464277	0.00965614531108591...
    0.00957887951808300	0.00950152956764568	0.00942651005452615...
    0.00935564857930749	0.00929009530525501	0.00923034252391390...
    0.00917624325277688	0.00912660693207717	0.00908001085044452...
    0.00903498466968650	0.00899025428193140	0.00894467938980312...
    0.00889742597289350	0.00884834714198970	0.00879751488542994...
    0.00874500476298096	0.00869110202828708	0.00863603331400137...
    0.00858005750933664	0.00852344082420545	0.00846584927026959...
    0.00840767340146492	0.00834942826152323	0.00829143343353489...
    0.00823425762595472	0.00817840148309757]';

%input wavelengths [nm]
lambda = (400:700)';

%output percentiles
PT = [10 90];

%partition adg spectra using partitioning model with UVVIS library
[lambda_out,adopt,agopt] = Partition_ADG_VIS(lambda, adg, PT);

%save inputs and outputs into an excel file
T1 = table(lambda,adg);
T2 = table(lambda_out,adopt(:,1),agopt(:,1),adopt(:,1)+agopt(:,1),adopt(:,2),agopt(:,2),adopt(:,2)+agopt(:,2),adopt(:,3),agopt(:,3),adopt(:,3)+agopt(:,3));
T1.Properties.VariableNames = {'Input Wavelength [nm]','Input adg [1/m]'};
T2.Properties.VariableNames = {'Output Wavelength [nm]','Output ad (Optimal) [1/m]','Output ag (Optimal) [1/m]','Output adg (Optimal) [1/m]','Output ad (10th) [1/m]','Output ag (10th) [1/m]','Output adg (10th) [1/m]','Output ad (90th) [1/m]','Output ag (90th) [1/m]','Output adg (90th) [1/m]'};
FormatOut = 'yyyymmdd';
outfile = ['Partition_ADG_VIS_Test_Run_' datestr(datetime,FormatOut)];
writetable(T1,outfile,'FileType','spreadsheet','Sheet','ADG_part_Input')
writetable(T2,outfile,'FileType','spreadsheet','Sheet','ADG_part_Output_VIS_UVext')
