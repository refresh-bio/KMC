/*
 *
 * nthash.hpp
 * Author: Hamid Mohamadi
 * Genome Sciences Centre,
 * British Columbia Cancer Agency
 */

#ifndef NT_HASH_H
#define NT_HASH_H

#include <limits>
#include <stdint.h>
#include <vector>

static const uint8_t RCconvertTab[256] = {
    255, 255, 255, 255, 255, 255, 255, 255, // 0..7
    255, 255, 255, 255, 255, 255, 255, 255, // 8..15
    255, 255, 255, 255, 255, 255, 255, 255, // 16..23
    255, 255, 255, 255, 255, 255, 255, 255, // 24..31
    255, 255, 255, 255, 255, 255, 255, 255, // 32..39
    255, 255, 255, 255, 255, 255, 255, 255, // 40..47
    255, 255, 255, 255, 255, 255, 255, 255, // 48..55
    255, 255, 255, 255, 255, 255, 255, 255, // 56..63
    255, 3, 255, 2, 255, 255, 255, 1, // 64..71
    255, 255, 255, 255, 255, 255, 255, 255, // 72..79
    255, 255, 255, 255, 0, 3, 255, 255, // 80..87
    255, 255, 255, 255, 255, 255, 255, 255, // 88..95
    255, 3, 255, 2, 255, 255, 255, 1, // 96..103
    255, 255, 255, 255, 255, 255, 255, 255, // 104..111
    255, 255, 255, 255, 0, 3, 255, 255, // 112..119
    255, 255, 255, 255, 255, 255, 255, 255, // 120..127
    255, 255, 255, 255, 255, 255, 255, 255, // 128..135
    255, 255, 255, 255, 255, 255, 255, 255, // 136..143
    255, 255, 255, 255, 255, 255, 255, 255, // 144..151
    255, 255, 255, 255, 255, 255, 255, 255, // 152..159
    255, 255, 255, 255, 255, 255, 255, 255, // 160..167
    255, 255, 255, 255, 255, 255, 255, 255, // 168..175
    255, 255, 255, 255, 255, 255, 255, 255, // 176..183
    255, 255, 255, 255, 255, 255, 255, 255, // 184..191
    255, 255, 255, 255, 255, 255, 255, 255, // 192..199
    255, 255, 255, 255, 255, 255, 255, 255, // 200..207
    255, 255, 255, 255, 255, 255, 255, 255, // 208..215
    255, 255, 255, 255, 255, 255, 255, 255, // 216..223
    255, 255, 255, 255, 255, 255, 255, 255, // 224..231
    255, 255, 255, 255, 255, 255, 255, 255, // 232..239
    255, 255, 255, 255, 255, 255, 255, 255, // 240..247
    255, 255, 255, 255, 255, 255, 255, 255  // 248..255
};

static const uint8_t convertTab[256] = {
    255, 255, 255, 255, 255, 255, 255, 255, // 0..7
    255, 255, 255, 255, 255, 255, 255, 255, // 8..15
    255, 255, 255, 255, 255, 255, 255, 255, // 16..23
    255, 255, 255, 255, 255, 255, 255, 255, // 24..31
    255, 255, 255, 255, 255, 255, 255, 255, // 32..39
    255, 255, 255, 255, 255, 255, 255, 255, // 40..47
    255, 255, 255, 255, 255, 255, 255, 255, // 48..55
    255, 255, 255, 255, 255, 255, 255, 255, // 56..63
    255, 0, 255, 1, 255, 255, 255, 2, // 64..71
    255, 255, 255, 255, 255, 255, 255, 255, // 72..79
    255, 255, 255, 255, 3, 0, 255, 255, // 80..87
    255, 255, 255, 255, 255, 255, 255, 255, // 88..95
    255, 0, 255, 1, 255, 255, 255, 2, // 96..103
    255, 255, 255, 255, 255, 255, 255, 255, // 104..111
    255, 255, 255, 255, 3, 0, 255, 255, // 112..119
    255, 255, 255, 255, 255, 255, 255, 255, // 120..127
    255, 255, 255, 255, 255, 255, 255, 255, // 128..135
    255, 255, 255, 255, 255, 255, 255, 255, // 136..143
    255, 255, 255, 255, 255, 255, 255, 255, // 144..151
    255, 255, 255, 255, 255, 255, 255, 255, // 152..159
    255, 255, 255, 255, 255, 255, 255, 255, // 160..167
    255, 255, 255, 255, 255, 255, 255, 255, // 168..175
    255, 255, 255, 255, 255, 255, 255, 255, // 176..183
    255, 255, 255, 255, 255, 255, 255, 255, // 184..191
    255, 255, 255, 255, 255, 255, 255, 255, // 192..199
    255, 255, 255, 255, 255, 255, 255, 255, // 200..207
    255, 255, 255, 255, 255, 255, 255, 255, // 208..215
    255, 255, 255, 255, 255, 255, 255, 255, // 216..223
    255, 255, 255, 255, 255, 255, 255, 255, // 224..231
    255, 255, 255, 255, 255, 255, 255, 255, // 232..239
    255, 255, 255, 255, 255, 255, 255, 255, // 240..247
    255, 255, 255, 255, 255, 255, 255, 255  // 248..255
};

static const uint64_t dimerTab[16] = {
5015898201438948509U, 5225361804584821669U, 6423762225589857229U, 5783394398799547583U,
6894017875502584557U, 5959461383092338133U, 4833978511655400893U, 5364573296520205007U,
9002561594443973180U, 8212239310050454788U, 6941810030513055084U, 7579897184553533982U,
7935738758488558809U, 7149836515649299425U, 8257540373175577481U, 8935100007508790523U
};

static const uint64_t trimerTab[64] = {
13237172352163388750U, 13451082378889146998U, 12324706752351386142U, 11704099346423635308U,
12503002411303846718U, 11573033083854154758U, 12770611021816489070U, 13284814289517544220U,
10286336837755622383U, 9500434588327378135U, 10554658215321236671U, 11177611689138066381U,
11245073286936829194U, 10454751004568891954U, 9274956656780491354U, 9930495270120774952U,
9498947889754972591U, 10289371588586147479U, 11487222103436658431U, 10812501148518244749U,
11088845979783725023U, 10735249574334615783U, 9609199230360475791U, 10105458452942995453U,
13447889238169808654U, 13238535845420384310U, 11968673763542288478U, 12645600078955589420U,
12136759312206930411U, 11922809957208297171U, 13031072242070652603U, 13668666814620918217U,
14219262150204358668U, 14433136993975185204U, 15703263506252408668U, 15026899868095529006U,
16097136083696541308U, 15167201938128040260U, 14113514427211577644U, 14608043031429815902U,
18169629015343943341U, 17383691583363408277U, 16185576633819064829U, 16859734366019948175U,
17215452794964541512U, 16425095330967072624U, 17460550829194815256U, 18101973914136232042U,
16197524846324948423U, 17136496960994620159U, 18190301010467109527U, 17660752969549176293U,
18084590689685816247U, 17861669045228104847U, 16591430392433501415U, 17233003275094786965U,
15689030113991676774U, 15321980360070757470U, 14196301091602199606U, 14727918144983470916U,
14660430141886012803U, 14297932370981794491U, 15550237822687034067U, 16044915679164358049U
};

static const uint64_t tetramerTab[256] = {
6047278271377325800U, 6842100033257738704U, 5716751207778949560U, 5058261232784932554U,
5322212292231585944U, 4955210659836481440U, 6153481158060361672U, 6630136099103187130U,
7683058811908681801U, 7460089081761259377U, 8513615477720831769U, 9169618076073996395U,
8669810821731892908U, 8451393064794886548U, 7271235746105367036U, 7894785163577458318U,
7461575445318369801U, 7680024275870068017U, 8878022265940976985U, 8237757801848291883U,
9060296013225843833U, 8116780716040188737U, 6991106539262573353U, 7521593563379047515U,
6845292839028968616U, 6045914992845185936U, 4775672622745250808U, 5413871935584767114U,
5490367161684853325U, 4695435745326017909U, 5803018666222232861U, 6480400171096490607U,
2381043025085637546U, 3175899973157948562U, 4445879008075678970U, 3807116472585741192U,
4268108881087626714U, 3901072061426881250U, 2847008385469766282U, 3379366782720458232U,
1763336001516006667U, 1540401457157816883U, 342666797974407771U, 983493939256405289U,
771890739233563630U, 553508169276984534U, 1589643033626739902U, 2263336780810576844U,
330722743541775969U, 688712796851212633U, 1742668713148160305U, 1245320973785726531U,
2208596672445898769U, 1422777727841816361U, 152919646732699457U, 826464124477841459U,
4460107693596700864U, 3530055095011467256U, 2403999925630162832U, 2899137386794791138U,
3398970977768160805U, 2464498338584432925U, 3716128830812494197U, 4248337413163712007U,
4264326372183459627U, 3906261395711551507U, 2851952150714671227U, 3383149429014333193U,
2386233046276708699U, 3172117876357805667U, 4441779805226941963U, 3801926588820052345U,
170684860043692426U, 1100671402695403186U, 2226926226858061530U, 1693589575942097320U,
1193606390847620975U, 2128144916583147607U, 876319371625685055U, 382305650241144653U,
1102545060664966090U, 168107437338776818U, 1437989166537956506U, 1915072878734195688U,
1548519783094789562U, 1757891215679916674U, 703889661060612842U, 46092416782165400U,
3908715595921208683U, 4262294307145226835U, 3064498623987880507U, 2585134797421409609U,
2661735585529691022U, 3019760716990469302U, 4055956603131813086U, 3543998858204232620U,
5317339067591416425U, 4959238909506745681U, 6157334207435046201U, 6635009461133220427U,
6051307208490845209U, 6837227221258447649U, 5711490920986878793U, 5054232433096901691U,
8122648135453742280U, 9052599496358476784U, 7782418148093113240U, 7307023562816214250U,
7095314801322056237U, 8029818144085865749U, 9137340041034366333U, 8622472983995947535U,
7806751516869674914U, 7011855109925922970U, 8137690373747335410U, 8757695200062998400U,
8531879593853721042U, 8898947385530005226U, 7700757522090507906U, 7186022138009770480U,
6135219772853324035U, 6358123720871388731U, 5304510851123850835U, 4682089562405882145U,
5182028715320330214U, 5400512630465816798U, 6580751683450298550U, 5923625422568720324U,
13124074928584983660U, 13491146941631638356U, 12293650504952193852U, 11816502978180760654U,
12399079312662682140U, 11604187204414436644U, 12730450818222161228U, 13388307479092468286U,
10327209524901530317U, 9388215691182564853U, 10657868830410829213U, 11137168911054473967U,
11357920004770333736U, 10414374197647485712U, 9306325182584103800U, 9818342344138146826U,
9386341947321596045U, 10329786896059045813U, 11455812913355464669U, 10924692575052363951U,
10984992149858150141U, 10766613702172592581U, 9568826821541020077U, 10208598699842184927U,
13488692655530571308U, 13126106942075820308U, 12072096584926548348U, 12605510244625659406U,
12249677498819492041U, 11882645355480553457U, 13062230760632229785U, 13556163143878539499U,
14178740190036597038U, 14545847390080448022U, 15599559227675164286U, 15067834145139579148U,
16065876409530435422U, 15270949115358734438U, 14000758968863088654U, 14640014089599289212U,
18281953465151117199U, 17342994818563569847U, 16217267316526477535U, 16746698532205467565U,
17255653680509032810U, 16312143059561297490U, 17564497017566543418U, 18061360711745100104U,
16237972021990524133U, 17023861349393640413U, 18293930539975648181U, 17619893477009409223U,
18115916316835994261U, 17757855915011241389U, 16704251839199542725U, 17200966263939144375U,
15576639675766950468U, 15362743113290245500U, 14164544455910714644U, 14841019967217601126U,
14620295210399335585U, 14410818688327658393U, 15446357621659116529U, 16085462927495578755U,
18237799192036655099U, 17294270664133710019U, 16258109964509321387U, 16773410497518403545U,
16657084189963477387U, 16875519862962278067U, 18127020052323321563U, 17507580374969491881U,
14153168177888129370U, 14515696771658964578U, 15624080140268688906U, 15110866744451150200U,
15466708232756051903U, 15833797605570023559U, 14563810316809509103U, 14085706539145691037U,
14517711175708869402U, 14150731501263563810U, 15402451490950456394U, 15899948742203982648U,
15224753927964908906U, 16019597712369578578U, 14983744703118572090U, 14310050713553640776U,
17296865610423782843U, 18235907873078829699U, 17055988043521714923U, 16561000163437350297U,
16340222631939670878U, 17283720110790814822U, 18338064546595415054U, 17805706452459078524U,
10375933128878629561U, 9432369415202180481U, 10612588863825479145U, 11105888166746317467U,
10794790039591648457U, 11013260899437695985U, 9905396050428550041U, 9228014311730625771U,
13154226096333843480U, 13516719503928509216U, 12264699899470662472U, 11768891770841246778U,
11836546934201131773U, 12203601119882644933U, 13328994472388527533U, 12798507759874630367U,
12277767672444305266U, 12068343612890878026U, 13176021535246260258U, 13816435502572994384U,
12705517425460601090U, 13640043170446921274U, 12460006250421962322U, 11929369723008524576U,
10597232027372843475U, 11387585128312430315U, 10351852510211364483U, 9713802769929286129U,
9357917249443839798U, 10143859113470169102U, 11342251114164164710U, 10664720106027613972U
};

// offset for the complement base in the random seeds table
const uint8_t cpOff = 0x07;

// shift for gerenerating multiple hash values
const int multiShift = 27;

// seed for gerenerating multiple hash values
static const uint64_t multiSeed = 0x90b45d39fb6da1fa;

// 64-bit random seeds corresponding to bases and their complements
static const uint64_t seedA = 0x3c8bfbb395c60474;
static const uint64_t seedC = 0x3193c18562a02b4c;
static const uint64_t seedG = 0x20323ed082572324;
static const uint64_t seedT = 0x295549f54be24456;
static const uint64_t seedN = 0x0000000000000000;

static const uint64_t seedTab[256] = {
    seedN, seedT, seedN, seedG, seedA, seedA, seedN, seedC, // 0..7
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 8..15
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 16..23
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 24..31
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 32..39
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 40..47
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 48..55
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 56..63
    seedN, seedA, seedN, seedC, seedN, seedN, seedN, seedG, // 64..71
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 72..79
    seedN, seedN, seedN, seedN, seedT, seedT, seedN, seedN, // 80..87
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 88..95
    seedN, seedA, seedN, seedC, seedN, seedN, seedN, seedG, // 96..103
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 104..111
    seedN, seedN, seedN, seedN, seedT, seedT, seedN, seedN, // 112..119
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 120..127
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 128..135
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 136..143
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 144..151
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 152..159
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 160..167
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 168..175
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 176..183
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 184..191
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 192..199
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 200..207
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 208..215
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 216..223
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 224..231
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 232..239
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 240..247
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN  // 248..255
};

static const uint64_t A33r[33] = {
    0x195c60474,0x12b8c08e9,0x571811d3,0xae3023a6,0x15c60474c,0xb8c08e99,0x171811d32,0xe3023a65,0x1c60474ca,0x18c08e995,0x11811d32b,0x3023a657,0x60474cae,0xc08e995c,0x1811d32b8,0x1023a6571,0x474cae3,0x8e995c6,0x11d32b8c,0x23a65718,0x474cae30,0x8e995c60,0x11d32b8c0,0x3a657181,0x74cae302,0xe995c604,0x1d32b8c08,0x1a6571811,0x14cae3023,0x995c6047,0x132b8c08e,0x6571811d,0xcae3023a
};

static const uint64_t A31l[31] = {
    0x3c8bfbb200000000,0x7917f76400000000,0xf22feec800000000,0xe45fdd9200000000,0xc8bfbb2600000000,0x917f764e00000000,0x22feec9e00000000,0x45fdd93c00000000,0x8bfbb27800000000,0x17f764f200000000,0x2feec9e400000000,0x5fdd93c800000000,0xbfbb279000000000,0x7f764f2200000000,0xfeec9e4400000000,0xfdd93c8a00000000,0xfbb2791600000000,0xf764f22e00000000,0xeec9e45e00000000,0xdd93c8be00000000,0xbb27917e00000000,0x764f22fe00000000,0xec9e45fc00000000,0xd93c8bfa00000000,0xb27917f600000000,0x64f22fee00000000,0xc9e45fdc00000000,0x93c8bfba00000000,0x27917f7600000000,0x4f22feec00000000,0x9e45fdd800000000
};

static const uint64_t C33r[33] = {
    0x162a02b4c,0xc5405699,0x18a80ad32,0x115015a65,0x2a02b4cb,0x54056996,0xa80ad32c,0x15015a658,0xa02b4cb1,0x140569962,0x80ad32c5,0x1015a658a,0x2b4cb15,0x569962a,0xad32c54,0x15a658a8,0x2b4cb150,0x569962a0,0xad32c540,0x15a658a80,0xb4cb1501,0x169962a02,0xd32c5405,0x1a658a80a,0x14cb15015,0x9962a02b,0x132c54056,0x658a80ad,0xcb15015a,0x1962a02b4,0x12c540569,0x58a80ad3,0xb15015a6
};

static const uint64_t C31l[31] = {
    0x3193c18400000000,0x6327830800000000,0xc64f061000000000,0x8c9e0c2200000000,0x193c184600000000,0x3278308c00000000,0x64f0611800000000,0xc9e0c23000000000,0x93c1846200000000,0x278308c600000000,0x4f06118c00000000,0x9e0c231800000000,0x3c18463200000000,0x78308c6400000000,0xf06118c800000000,0xe0c2319200000000,0xc184632600000000,0x8308c64e00000000,0x6118c9e00000000,0xc23193c00000000,0x1846327800000000,0x308c64f000000000,0x6118c9e000000000,0xc23193c000000000,0x8463278200000000,0x8c64f0600000000,0x118c9e0c00000000,0x23193c1800000000,0x4632783000000000,0x8c64f06000000000,0x18c9e0c200000000
};


static const uint64_t G33r[33] = {
    0x82572324,0x104ae4648,0x95c8c91,0x12b91922,0x25723244,0x4ae46488,0x95c8c910,0x12b919220,0x57232441,0xae464882,0x15c8c9104,0xb9192209,0x172324412,0xe4648825,0x1c8c9104a,0x191922095,0x12324412b,0x46488257,0x8c9104ae,0x11922095c,0x324412b9,0x64882572,0xc9104ae4,0x1922095c8,0x124412b91,0x48825723,0x9104ae46,0x122095c8c,0x4412b919,0x88257232,0x1104ae464,0x2095c8c9,0x412b9192
};

static const uint64_t G31l[31] = {
    0x20323ed000000000,0x40647da000000000,0x80c8fb4000000000,0x191f68200000000,0x323ed0400000000,0x647da0800000000,0xc8fb41000000000,0x191f682000000000,0x323ed04000000000,0x647da08000000000,0xc8fb410000000000,0x91f6820200000000,0x23ed040600000000,0x47da080c00000000,0x8fb4101800000000,0x1f68203200000000,0x3ed0406400000000,0x7da080c800000000,0xfb41019000000000,0xf682032200000000,0xed04064600000000,0xda080c8e00000000,0xb410191e00000000,0x6820323e00000000,0xd040647c00000000,0xa080c8fa00000000,0x410191f600000000,0x820323ec00000000,0x40647da00000000,0x80c8fb400000000,0x10191f6800000000
};

static const uint64_t T33r[33] = {
    0x14be24456,0x97c488ad,0x12f89115a,0x5f1222b5,0xbe24456a,0x17c488ad4,0xf89115a9,0x1f1222b52,0x1e24456a5,0x1c488ad4b,0x189115a97,0x11222b52f,0x24456a5f,0x488ad4be,0x9115a97c,0x1222b52f8,0x4456a5f1,0x88ad4be2,0x1115a97c4,0x22b52f89,0x456a5f12,0x8ad4be24,0x115a97c48,0x2b52f891,0x56a5f122,0xad4be244,0x15a97c488,0xb52f8911,0x16a5f1222,0xd4be2445,0x1a97c488a,0x152f89115,0xa5f1222b
};

static const uint64_t T31l[31] = {
    0x295549f400000000,0x52aa93e800000000,0xa55527d000000000,0x4aaa4fa200000000,0x95549f4400000000,0x2aa93e8a00000000,0x55527d1400000000,0xaaa4fa2800000000,0x5549f45200000000,0xaa93e8a400000000,0x5527d14a00000000,0xaa4fa29400000000,0x549f452a00000000,0xa93e8a5400000000,0x527d14aa00000000,0xa4fa295400000000,0x49f452aa00000000,0x93e8a55400000000,0x27d14aaa00000000,0x4fa2955400000000,0x9f452aa800000000,0x3e8a555200000000,0x7d14aaa400000000,0xfa29554800000000,0xf452aa9200000000,0xe8a5552600000000,0xd14aaa4e00000000,0xa295549e00000000,0x452aa93e00000000,0x8a55527c00000000,0x14aaa4fa00000000
};

static const uint64_t N33r[33] = {
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN
};

static const uint64_t N31l[31] = {
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN
};


static const uint64_t *msTab33r[256] = {
    N33r, T33r, N33r, G33r, A33r, A33r, N33r, C33r, // 0..7
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 8..15
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 16..23
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 24..31
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 32..39
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 40..47
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 48..55
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 56..63
    N33r, A33r, N33r, C33r, N33r, N33r, N33r, G33r, // 64..71
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 72..79
    N33r, N33r, N33r, N33r, T33r, T33r, N33r, N33r, // 80..87
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 88..95
    N33r, A33r, N33r, C33r, N33r, N33r, N33r, G33r, // 96..103
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 104..111
    N33r, N33r, N33r, N33r, T33r, T33r, N33r, N33r, // 112..119
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 120..127
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 128..135
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 136..143
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 144..151
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 152..159
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 160..167
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 168..175
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 176..183
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 184..191
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 192..199
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 200..207
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 208..215
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 216..223
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 224..231
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 232..239
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 240..247
    N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r  // 248..255
};

static const uint64_t *msTab31l[256] = {
    N31l, T31l, N31l, G31l, A31l, A31l, N31l, C31l, // 0..7
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 8..15
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 16..23
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 24..31
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 32..39
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 40..47
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 48..55
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 56..63
    N31l, A31l, N31l, C31l, N31l, N31l, N31l, G31l, // 64..71
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 72..79
    N31l, N31l, N31l, N31l, T31l, T31l, N31l, N31l, // 80..87
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 88..95
    N31l, A31l, N31l, C31l, N31l, N31l, N31l, G31l, // 96..103
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 104..111
    N31l, N31l, N31l, N31l, T31l, T31l, N31l, N31l, // 112..119
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 120..127
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 128..135
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 136..143
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 144..151
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 152..159
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 160..167
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 168..175
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 176..183
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 184..191
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 192..199
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 200..207
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 208..215
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 216..223
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 224..231
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 232..239
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 240..247
    N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l  // 248..255
};

// rotate "v" to the left 1 position
inline uint64_t rol1(const uint64_t v) {
    return (v << 1) | (v >> 63);
}

// rotate "v" to the left x position
inline uint64_t rolx(const uint64_t v, const unsigned x) {
    return (v << x) | (v >> (64 - x));
}

// rotate "v" to the right by 1 position
inline uint64_t ror1(const uint64_t v) {
    return (v >> 1) | (v << 63);
}

// rotate 31-left bits of "v" to the left by "s" positions
inline uint64_t rol31(const uint64_t v, unsigned s) {
    s%=31;
    return ((v << s) | (v >> (31 - s))) & 0x7FFFFFFF;
}

// rotate 33-right bits of "v" to the left by "s" positions
inline uint64_t rol33(const uint64_t v, unsigned s) {
    s%=33;
    return ((v << s) | (v >> (33 - s))) & 0x1FFFFFFFF;
}

// swap bit 0 with bit 33 in "v"
inline uint64_t swapbits033(const uint64_t v) {
    uint64_t x = (v ^ (v >> 33)) & 1;
    return v ^ (x | (x << 33));
}

// swap bit 32 with bit 63 in "v"
inline uint64_t swapbits3263(const uint64_t v) {
    uint64_t x = ((v >> 32) ^ (v >> 63)) & 1;
    return v ^ ((x << 32) | (x << 63));
}

inline uint64_t swapxbits033(const uint64_t v, const unsigned x) {
    uint64_t y = (v ^ (v >> 33)) & ((std::numeric_limits<uint64_t>::max)() >> (64 - x));
    return v ^ (y | (y << 33));
}

// forward-strand hash value of the base kmer, i.e. fhval(kmer_0)
inline uint64_t NTF64(const char * kmerSeq, const unsigned k) {
    uint64_t hVal=0;
    for (unsigned i=0; i<k/4; i++) {
        hVal = rolx(hVal, 4);
        hVal = swapxbits033(hVal, 4);
		uint8_t currOffSet = 4 * i;
		uint8_t tetramerLoc = 64 * convertTab[(unsigned char)kmerSeq[currOffSet]] + 16 * convertTab[(unsigned char)kmerSeq[currOffSet + 1]] + 4 * convertTab[(unsigned char)kmerSeq[currOffSet + 2]] +  convertTab[(unsigned char)kmerSeq[currOffSet + 3]];
        hVal ^= tetramerTab[tetramerLoc];
    }
	unsigned remainder = k % 4;
    hVal = rolx(hVal, remainder);
    hVal = swapxbits033(hVal, remainder);
	if (remainder == 3) {
		uint8_t trimerLoc = 16 * convertTab[(unsigned char)kmerSeq[k - 3]] + 4 * convertTab[(unsigned char)kmerSeq[k - 2]] +  convertTab[(unsigned char)kmerSeq[k - 1]];
        hVal ^= trimerTab[trimerLoc];
	} else if (remainder == 2) {
		uint8_t dimerLoc = 4 * convertTab[(unsigned char)kmerSeq[k - 2]] +  convertTab[(unsigned char)kmerSeq[k - 1]];
        hVal ^= dimerTab[dimerLoc];
	}  else if (remainder == 1) {
		hVal ^= seedTab[(unsigned char)kmerSeq[k - 1]];
	}
    return hVal;
}

// reverse-strand hash value of the base kmer, i.e. rhval(kmer_0)
inline uint64_t NTR64(const char * kmerSeq, const unsigned k) {
    uint64_t hVal=0;
	unsigned remainder = k % 4;
	if (remainder == 3) {
		uint8_t trimerLoc = 16 * RCconvertTab[(unsigned char)kmerSeq[k - 1]] + 4 * RCconvertTab[(unsigned char)kmerSeq[k - 2]] +  RCconvertTab[(unsigned char)kmerSeq[k - 3]];
        hVal ^= trimerTab[trimerLoc];
	} else if (remainder == 2) {
		uint8_t dimerLoc = 4 * RCconvertTab[(unsigned char)kmerSeq[k - 1]] +  RCconvertTab[(unsigned char)kmerSeq[k - 2]];
        hVal ^= dimerTab[dimerLoc];
	} else if (remainder == 1) {
		hVal ^= seedTab[(unsigned char)kmerSeq[k-1]&cpOff];
	}
    for(unsigned i=0; i<k/4; i++) {
        hVal = rolx(hVal, 4);
        hVal = swapxbits033(hVal, 4);
		uint8_t currOffSet = 4 * (k / 4 - i) - 1;
		uint8_t tetramerLoc= 64 * RCconvertTab[(unsigned char)kmerSeq[currOffSet]] + 16 * RCconvertTab[(unsigned char)kmerSeq[currOffSet - 1]] + 4 * RCconvertTab[(unsigned char)kmerSeq[currOffSet - 2]] +  RCconvertTab[(unsigned char)kmerSeq[currOffSet - 3]];
		hVal ^= tetramerTab[tetramerLoc];
    }
    return hVal;
}

// forward-strand ntHash for sliding k-mers
inline uint64_t NTF64(const uint64_t fhVal, const unsigned k, const unsigned char charOut, const unsigned char charIn) {
    uint64_t hVal = rol1(fhVal);
    hVal = swapbits033(hVal);
    hVal ^= seedTab[charIn];
    hVal ^= (msTab31l[charOut][k%31] | msTab33r[charOut][k%33]);
    return hVal;
}

// reverse-complement ntHash for sliding k-mers
inline uint64_t NTR64(const uint64_t rhVal, const unsigned k, const unsigned char charOut, const unsigned char charIn) {
    uint64_t hVal = rhVal ^ (msTab31l[charIn&cpOff][k%31] | msTab33r[charIn&cpOff][k%33]);
    hVal ^= seedTab[charOut&cpOff];
    hVal = ror1(hVal);
    hVal = swapbits3263(hVal);
    return hVal;
}

// canonical ntBase
inline uint64_t NTC64(const char * kmerSeq, const unsigned k) {
    uint64_t fhVal=0, rhVal=0;
    fhVal=NTF64(kmerSeq, k);
    rhVal=NTR64(kmerSeq, k);
    return (rhVal<fhVal)? rhVal : fhVal;
}

// canonical ntHash
inline uint64_t NTC64(const char * kmerSeq, const unsigned k, uint64_t& fhVal, uint64_t& rhVal) {
    fhVal = NTF64(kmerSeq, k);
    rhVal = NTR64(kmerSeq, k);
    return (rhVal<fhVal)? rhVal : fhVal;
}

// canonical ntHash for sliding k-mers
inline uint64_t NTC64(const unsigned char charOut, const unsigned char charIn, const unsigned k, uint64_t& fhVal, uint64_t& rhVal) {
    fhVal = NTF64(fhVal, k, charOut, charIn);
    rhVal = NTR64(rhVal, k, charOut, charIn);
    return (rhVal<fhVal)? rhVal : fhVal;
}

// forward-strand ntHash for sliding k-mers to the left
inline uint64_t NTF64L(const uint64_t rhVal, const unsigned k, const unsigned char charOut, const unsigned char charIn) {
    uint64_t hVal = rhVal ^ (msTab31l[charIn][k%31] | msTab33r[charIn][k%33]);
    hVal ^= seedTab[charOut];
    hVal = ror1(hVal);
    hVal = swapbits3263(hVal);
    return hVal;
}

// reverse-complement ntHash for sliding k-mers to the left
inline uint64_t NTR64L(const uint64_t fhVal, const unsigned k, const unsigned char charOut, const unsigned char charIn) {
    uint64_t hVal = rol1(fhVal);
    hVal = swapbits033(hVal);
    hVal ^= seedTab[charIn&cpOff];
    hVal ^= (msTab31l[charOut&cpOff][k%31] | msTab33r[charOut&cpOff][k%33]);
    return hVal;
}

// canonical ntHash for sliding k-mers to the left
inline uint64_t NTC64L(const unsigned char charOut, const unsigned char charIn, const unsigned k, uint64_t& fhVal, uint64_t& rhVal) {
    fhVal = NTF64L(fhVal, k, charOut, charIn);
    rhVal = NTR64L(rhVal, k, charOut, charIn);
    return (rhVal<fhVal)? rhVal : fhVal;
}

// ntBase with seeding option
inline uint64_t NTF64(const char * kmerSeq, const unsigned k, const unsigned seed) {
    uint64_t hVal=NTF64(kmerSeq, k);
    if(seed==0) return hVal;
    hVal *= seed ^ k * multiSeed;
    hVal ^= hVal >> multiShift;
    return hVal;
}

// canonical ntBase with seeding option
inline uint64_t NTC64(const char * kmerSeq, const unsigned k, const unsigned seed) {
    uint64_t hVal = NTC64(kmerSeq,k);
    if(seed==0) return hVal;
    hVal *= seed ^ k * multiSeed;
    hVal ^= hVal >> multiShift;
    return hVal;
}

// multihash ntHash, ntBase
inline void NTM64(const char * kmerSeq, const unsigned k, const unsigned m, uint64_t *hVal) {
    uint64_t bVal=0, tVal=0;
    bVal = NTF64(kmerSeq, k);
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

// one extra hash for given base hash
inline uint64_t NTE64(const uint64_t hVal, const unsigned k, const unsigned i) {
    uint64_t tVal = hVal;
    tVal *= (i ^ k * multiSeed);
    tVal ^= tVal >> multiShift;
    return tVal;
}

// multihash ntHash for sliding k-mers
inline void NTM64(const unsigned char charOut, const unsigned char charIn, const unsigned k, const unsigned m, uint64_t *hVal) {
    uint64_t bVal=0, tVal=0;
    bVal = NTF64(hVal[0], k, charOut, charIn);
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

// canonical multihash ntBase
inline void NTMC64(const char * kmerSeq, const unsigned k, const unsigned m, uint64_t *hVal) {
    uint64_t bVal=0, tVal=0;
    bVal = NTC64(kmerSeq, k);
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

// canonical multihash ntHash
inline void NTMC64(const char * kmerSeq, const unsigned k, const unsigned m, uint64_t& fhVal, uint64_t& rhVal, uint64_t *hVal) {
    uint64_t bVal=0, tVal=0;
    bVal = NTC64(kmerSeq, k, fhVal, rhVal);
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

// canonical multihash ntHash for sliding k-mers
inline void NTMC64(const unsigned char charOut, const unsigned char charIn, const unsigned k, const unsigned m, uint64_t& fhVal, uint64_t& rhVal, uint64_t *hVal) {
    uint64_t bVal=0, tVal=0;
    bVal = NTC64(charOut, charIn, k, fhVal, rhVal);
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

/*
 * ignoring k-mers containing nonACGT using ntHash function
*/

// canonical ntBase
inline bool NTC64(const char *kmerSeq, const unsigned k, uint64_t& hVal, unsigned& locN) {
    hVal=0;
    locN=0;
    uint64_t fhVal=0,rhVal=0;
    for(int i=k-1; i>=0; i--) {
        if(seedTab[(unsigned char)kmerSeq[i]]==seedN) {
            locN=i;
            return false;
        }
        fhVal = rol1(fhVal);
        fhVal = swapbits033(fhVal);
        fhVal ^= seedTab[(unsigned char)kmerSeq[k-1-i]];

        rhVal = rol1(rhVal);
        rhVal = swapbits033(rhVal);
        rhVal ^= seedTab[(unsigned char)kmerSeq[i]&cpOff];
    }
    hVal = (rhVal<fhVal)? rhVal : fhVal;
    return true;
}

// canonical multihash ntBase
inline bool NTMC64(const char *kmerSeq, const unsigned k, const unsigned m, unsigned& locN, uint64_t* hVal) {
    uint64_t bVal=0, tVal=0, fhVal=0, rhVal=0;
    locN=0;
    for(int i=k-1; i>=0; i--) {
        if(seedTab[(unsigned char)kmerSeq[i]]==seedN) {
            locN=i;
            return false;
        }
        fhVal = rol1(fhVal);
        fhVal = swapbits033(fhVal);
        fhVal ^= seedTab[(unsigned char)kmerSeq[k-1-i]];

        rhVal = rol1(rhVal);
        rhVal = swapbits033(rhVal);
        rhVal ^= seedTab[(unsigned char)kmerSeq[i]&cpOff];
    }
    bVal = (rhVal<fhVal)? rhVal : fhVal;
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
    return true;
}

// canonical ntHash
inline bool NTC64(const char *kmerSeq, const unsigned k, uint64_t& fhVal, uint64_t& rhVal, uint64_t& hVal, unsigned& locN) {
    hVal=fhVal=rhVal=0;
    locN=0;
    for(int i=k-1; i>=0; i--) {
        if(seedTab[(unsigned char)kmerSeq[i]]==seedN) {
            locN=i;
            return false;
        }
        fhVal = rol1(fhVal);
        fhVal = swapbits033(fhVal);
        fhVal ^= seedTab[(unsigned char)kmerSeq[k-1-i]];

        rhVal = rol1(rhVal);
        rhVal = swapbits033(rhVal);
        rhVal ^= seedTab[(unsigned char)kmerSeq[i]&cpOff];
    }
    hVal = (rhVal<fhVal)? rhVal : fhVal;
    return true;
}

// canonical multihash ntHash
inline bool NTMC64(const char *kmerSeq, const unsigned k, const unsigned m, uint64_t& fhVal, uint64_t& rhVal, unsigned& locN, uint64_t* hVal) {
    fhVal=rhVal=0;
    uint64_t bVal=0, tVal=0;
    locN=0;
    for(int i=k-1; i>=0; i--) {
        if(seedTab[(unsigned char)kmerSeq[i]]==seedN) {
            locN=i;
            return false;
        }
        fhVal = rol1(fhVal);
        fhVal = swapbits033(fhVal);
        fhVal ^= seedTab[(unsigned char)kmerSeq[k-1-i]];

        rhVal = rol1(rhVal);
        rhVal = swapbits033(rhVal);
        rhVal ^= seedTab[(unsigned char)kmerSeq[i]&cpOff];
    }
    bVal = (rhVal<fhVal)? rhVal : fhVal;
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
    return true;
}

// strand-aware canonical multihash ntHash
inline bool NTMC64(const char *kmerSeq, const unsigned k, const unsigned m, uint64_t& fhVal, uint64_t& rhVal, unsigned& locN, uint64_t* hVal, bool& hStn) {
    fhVal=rhVal=0;
    uint64_t bVal=0, tVal=0;
    locN=0;
    for(int i=k-1; i>=0; i--) {
        if(seedTab[(unsigned char)kmerSeq[i]]==seedN) {
            locN=i;
            return false;
        }
        fhVal = rol1(fhVal);
        fhVal = swapbits033(fhVal);
        fhVal ^= seedTab[(unsigned char)kmerSeq[k-1-i]];
        
        rhVal = rol1(rhVal);
        rhVal = swapbits033(rhVal);
        rhVal ^= seedTab[(unsigned char)kmerSeq[i]&cpOff];
    }
    hStn = rhVal<fhVal;
    bVal = hStn? rhVal : fhVal;
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
    return true;
}

// starnd-aware canonical multihash ntHash for sliding k-mers
inline void NTMC64(const unsigned char charOut, const unsigned char charIn, const unsigned k, const unsigned m, uint64_t& fhVal, uint64_t& rhVal, uint64_t *hVal, bool &hStn) {
    uint64_t bVal=0, tVal=0;
    bVal = NTC64(charOut, charIn, k, fhVal, rhVal);
    hStn = rhVal<fhVal;
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

// masking canonical ntHash using spaced seed pattern
inline uint64_t maskHash(uint64_t &fkVal, uint64_t &rkVal, const char * seedSeq, const char * kmerSeq, const unsigned k) {
    uint64_t fsVal=fkVal, rsVal=rkVal;
    for(unsigned i=0; i<k; i++) {
        if(seedSeq[i]!='1') {
            fsVal ^= (msTab31l[(unsigned char)kmerSeq[i]][(k-1-i)%31] | msTab33r[(unsigned char)kmerSeq[i]][(k-1-i)%33]);
            rsVal ^= (msTab31l[(unsigned char)kmerSeq[i]&cpOff][i%31] | msTab33r[(unsigned char)kmerSeq[i]&cpOff][i%33]);
            
        }
    }
    return (rsVal<fsVal)? rsVal : fsVal;
}

// spaced seed ntHash for base kmer, i.e. fhval(kmer_0)
inline uint64_t NTS64(const char * kmerSeq, const std::vector<bool> &seed, const unsigned k, uint64_t &hVal) {
    hVal=0;
    uint64_t sVal = 0;
    for(unsigned i=0; i<k; i++) {
        hVal = rol1(hVal);
        hVal = swapbits033(hVal);
        sVal = hVal;
        hVal ^= seedTab[(unsigned char)kmerSeq[i]];
        if(seed[i]) sVal = hVal;
    }
    return sVal;
}

// spaced seed ntHash for sliding k-mers
inline uint64_t NTS64(const char * kmerSeq, const std::vector<bool> &seed, const unsigned char charOut, const unsigned char charIn, const unsigned k, uint64_t &hVal) {
    hVal = NTF64(hVal,k,charOut,charIn);
    uint64_t sVal = hVal;
    for(unsigned i=0; i<k; i++)
        if(!seed[i]) {
            sVal ^= (msTab31l[(unsigned char)kmerSeq[i]][k%31] | msTab33r[(unsigned char)kmerSeq[i]][k%33]);
        }
    return sVal;
}

// strand-aware multihash spaced seed ntHash
inline bool NTMS64(const char *kmerSeq, const std::vector<std::vector<unsigned> > &seedSeq, const unsigned k, const unsigned m, uint64_t& fhVal, uint64_t& rhVal, unsigned& locN, uint64_t* hVal, bool *hStn) {
    fhVal=rhVal=0;
    locN=0;
    for(int i=k-1; i>=0; i--) {
        if(seedTab[(unsigned char)kmerSeq[i]]==seedN) {
            locN=i;
            return false;
        }
        fhVal = rol1(fhVal);
        fhVal = swapbits033(fhVal);
        fhVal ^= seedTab[(unsigned char)kmerSeq[k-1-i]];
        
        rhVal = rol1(rhVal);
        rhVal = swapbits033(rhVal);
        rhVal ^= seedTab[(unsigned char)kmerSeq[i]&cpOff];
    }
    
    for(unsigned j=0; j<m; j++) {
        uint64_t fsVal=fhVal, rsVal=rhVal;
        for(std::vector<unsigned>::const_iterator i=seedSeq[j].begin(); i!=seedSeq[j].end(); ++i) {
            fsVal ^= (msTab31l[(unsigned char)kmerSeq[*i]][(k-1-*i)%31] | msTab33r[(unsigned char)kmerSeq[*i]][(k-1-*i)%33]);
            rsVal ^= (msTab31l[(unsigned char)kmerSeq[*i]&cpOff][*i%31] | msTab33r[(unsigned char)kmerSeq[*i]&cpOff][*i%33]);
        }
        hStn[j] = rsVal<fsVal;
        hVal[j] = hStn[j]? rsVal : fsVal;
    }
    return true;
}

// strand-aware multihash spaced seed ntHash for sliding k-mers
inline void NTMS64(const char *kmerSeq, const std::vector<std::vector<unsigned> > &seedSeq, const unsigned char charOut, const unsigned char charIn, const unsigned k, const unsigned m, uint64_t& fhVal, uint64_t& rhVal, uint64_t *hVal, bool *hStn) {
    fhVal = NTF64(fhVal,k,charOut,charIn);
    rhVal = NTR64(rhVal,k,charOut,charIn);
    for(unsigned j=0; j<m; j++) {
        uint64_t fsVal=fhVal, rsVal=rhVal;
        for(std::vector<unsigned>::const_iterator i=seedSeq[j].begin(); i!=seedSeq[j].end(); ++i) {
            fsVal ^= (msTab31l[(unsigned char)kmerSeq[*i]][(k-1-*i)%31] | msTab33r[(unsigned char)kmerSeq[*i]][(k-1-*i)%33]);
            rsVal ^= (msTab31l[(unsigned char)kmerSeq[*i]&cpOff][*i%31] | msTab33r[(unsigned char)kmerSeq[*i]&cpOff][*i%33]);;
        }
        hStn[j] = rsVal<fsVal;
        hVal[j] = hStn[j]? rsVal : fsVal;
    }
}

// strand-aware multihash spaced seed ntHash with multiple hashes per seed
inline bool NTMSM64(const char *kmerSeq, const std::vector<std::vector<unsigned> > &seedSeq, const unsigned k, const unsigned m, const unsigned m2,
    uint64_t& fhVal, uint64_t& rhVal, unsigned& locN, uint64_t* hVal, bool *hStn)
{
    fhVal=rhVal=0;
    locN=0;
    for(int i=k-1; i>=0; i--) {
        if(seedTab[(unsigned char)kmerSeq[i]]==seedN) {
            locN=i;
            return false;
        }
        fhVal = rol1(fhVal);
        fhVal = swapbits033(fhVal);
        fhVal ^= seedTab[(unsigned char)kmerSeq[k-1-i]];

        rhVal = rol1(rhVal);
        rhVal = swapbits033(rhVal);
        rhVal ^= seedTab[(unsigned char)kmerSeq[i]&cpOff];
    }

    for(unsigned j=0; j<m; j++) {
        uint64_t fsVal=fhVal, rsVal=rhVal;
        for(std::vector<unsigned>::const_iterator i=seedSeq[j].begin(); i!=seedSeq[j].end(); ++i) {
            fsVal ^= (msTab31l[(unsigned char)kmerSeq[*i]][(k-1-*i)%31] | msTab33r[(unsigned char)kmerSeq[*i]][(k-1-*i)%33]);
            rsVal ^= (msTab31l[(unsigned char)kmerSeq[*i]&cpOff][*i%31] | msTab33r[(unsigned char)kmerSeq[*i]&cpOff][*i%33]);
        }
        hStn[j * m2] = rsVal<fsVal;
        hVal[j * m2] = hStn[j * m2]? rsVal : fsVal;
        for(unsigned j2=1; j2<m2; j2++) {
            uint64_t tVal = hVal[j * m2] * (j2 ^ k * multiSeed);
            tVal ^= tVal >> multiShift;
            hStn[j * m2 + j2] = hStn[j * m2];
            hVal[j * m2 + j2] = tVal;
        }
    }
    return true;
}

// strand-aware multihash spaced seed ntHash for sliding k-mers with multiple hashes per seed
inline void NTMSM64(const char *kmerSeq, const std::vector<std::vector<unsigned> > &seedSeq, const unsigned char charOut, const unsigned char charIn,
const unsigned k, const unsigned m, const unsigned m2, uint64_t& fhVal, uint64_t& rhVal, uint64_t *hVal, bool *hStn)
{
    fhVal = NTF64(fhVal,k,charOut,charIn);
    rhVal = NTR64(rhVal,k,charOut,charIn);
    for(unsigned j=0; j<m; j++) {
        uint64_t fsVal=fhVal, rsVal=rhVal;
        for(std::vector<unsigned>::const_iterator i=seedSeq[j].begin(); i!=seedSeq[j].end(); ++i) {
            fsVal ^= (msTab31l[(unsigned char)kmerSeq[*i]][(k-1-*i)%31] | msTab33r[(unsigned char)kmerSeq[*i]][(k-1-*i)%33]);
            rsVal ^= (msTab31l[(unsigned char)kmerSeq[*i]&cpOff][*i%31] | msTab33r[(unsigned char)kmerSeq[*i]&cpOff][*i%33]);
        }
        hStn[j * m2] = rsVal<fsVal;
        hVal[j * m2] = hStn[j * m2]? rsVal : fsVal;
        for(unsigned j2=1; j2<m2; j2++) {
            uint64_t tVal = hVal[j * m2] * (j2 ^ k * multiSeed);
            tVal ^= tVal >> multiShift;
            hStn[j * m2 + j2] = hStn[j * m2];
            hVal[j * m2 + j2] = tVal;
        }
    }
}


#endif