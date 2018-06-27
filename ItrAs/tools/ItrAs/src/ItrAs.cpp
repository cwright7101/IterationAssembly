#include <ItrAs.hpp>
#include "karect/karect.h"
using namespace std;

/********************************************************************************/

// We define some constant strings for names of command line parameters
static const char* STR_MINK             = "-mink"             ;
static const char* STR_MAXK             = "-maxk"             ;
static const char* STR_STEP             = "-step"             ;
static const char* STR_MINCONTIG        = "-min_contig"       ;

static const char* STR_TRAVERSAL_KIND   = "-traversal"        ;
static const char* STR_FASTA_LINE_SIZE  = "-fasta-line"       ;
static const char* STR_KEEP_ISOLATED    = "-keep-isolated"    ;

static const char* progressFormat0      = "ItrAs : assembly"  ;

static const char* STR_BOOL_PRECORRECT  = "-pre_correct"      ;
static const char* STR_CELLTYPE         = "-celltype"         ;
static const char* STR_MATCHTYPE        = "-matchtype"        ;
static const char* STR_AGGRESIVENESS    = "-aggressive"       ;
static const char* STR_PRECOR_STAGES    = "-numstages"        ;
static const char* STR_PRECOR_TRIM      = "-trim_reads"       ;

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
ItrAs::ItrAs ()  : Tool ("ItrAs")
{
  // We add some custom arguments for command line interface
  setParser (new OptionsParser ("ItrAs")); 
  getParser()->push_front (new OptionOneParam (STR_URI_INPUT,"input reads (fasta/fastq/compressed) or hdf5 file",   false));
  getParser()->push_back(new OptionOneParam("-out", "output directory", false, "out"));

  OptionsParser* precorrectionParser = new OptionsParser ("pre-correction");
  precorrectionParser->push_front(new OptionOneParam(STR_BOOL_PRECORRECT, "Bool used to turn on(1)/off(0) read precorrection, default is off", false, "0"));
  precorrectionParser->push_front(new OptionOneParam(STR_CELLTYPE, "[haploid|diploid] : Specify the cell type.(default diploid)", false, "diploid"));
  precorrectionParser->push_front(new OptionOneParam(STR_MATCHTYPE, "[edit|hamming|insdel] : Specify the matching type. 'hamming' allows substitution errors only. 'edit' allows insertions, deletion, and substitutions with equal costs. 'insdel' is the same as 'edit', but the cost of substitutions is doubled. Use 'hamming' for Illumina datasets, and 'edit' for 454 datasets.", false, "hamming"));
  precorrectionParser->push_front(new OptionOneParam(STR_AGGRESIVENESS, "Specify the aggressiveness towards error correction", false, "0.42"));
  precorrectionParser->push_front(new OptionOneParam(STR_PRECOR_STAGES, "Specify the number of stages (1 or 2) for read precorrection", false, "1"));
  precorrectionParser->push_front(new OptionOneParam(STR_PRECOR_TRIM, "Specify whether or not to trim corrected reads. If assembler needs fixed length reads, set to no.", false, "no"));
  getParser()->push_back (precorrectionParser);

  OptionsParser* assemblyParser = new OptionsParser ("assembly");
  assemblyParser->push_front (new OptionOneParam (STR_FASTA_LINE_SIZE, "number of nucleotides per line in fasta output (0 means one line)",  false, "0"));
  assemblyParser->push_front (new OptionOneParam (STR_TRAVERSAL_KIND,  "traversal type ('contig', 'unitig')", false,  "contig"  ));
  assemblyParser->push_front (new OptionNoParam  (STR_KEEP_ISOLATED,   "keep short (<= max(2k, 150 bp)) isolated output sequences", false));
  assemblyParser->push_front (new OptionOneParam(STR_MINK, "min k value", false, "25"));
  assemblyParser->push_front (new OptionOneParam(STR_MAXK, "max k value", false, "95"));
  assemblyParser->push_front (new OptionOneParam(STR_STEP, "step value for iteration of k values", false, "1"));
  assemblyParser->push_front (new OptionOneParam(STR_MINCONTIG, "min contig length", false, "100"));
  getParser()->push_back (assemblyParser);

  OptionsParser* simplificationsParser = new OptionsParser ("graph simplifications");
  Simplifications<GraphUnitigsTemplate<32>,NodeGU,EdgeGU> graphSimplifications(nullptr, 1, false);
  simplificationsParser->push_back (new OptionNoParam  ("-no-bulge-removal", "ask to not perform bulge removal", false));
  simplificationsParser->push_back (new OptionNoParam  ("-no-tip-removal",   "ask to not perform tip removal", false));
  simplificationsParser->push_back (new OptionNoParam  ("-no-ec-removal",   "ask to not perform erroneous connection removal", false));
  simplificationsParser->push_back (new OptionOneParam ("-tip-len-topo-kmult", "remove all tips of length <= k * X bp",  false, to_string(graphSimplifications._tipLen_Topo_kMult)));
  simplificationsParser->push_back (new OptionOneParam ("-tip-len-rctc-kmult", "remove tips that pass coverage criteria, of length <= k * X bp",  false, to_string(graphSimplifications._tipLen_RCTC_kMult)));
  simplificationsParser->push_back (new OptionOneParam ("-tip-rctc-cutoff",    "tip relative coverage coefficient: mean coverage of neighbors >  X * tip coverage",  false, to_string(graphSimplifications._tipRCTCcutoff)));
  simplificationsParser->push_back (new OptionOneParam ("-bulge-len-kmult",    "bulges shorter than k*X bp are candidate to be removed",  false, to_string(graphSimplifications._bulgeLen_kMult)));
  simplificationsParser->push_back (new OptionOneParam ("-bulge-len-kadd",     "bulges shorter than k+X bp are candidate to be removed",  false, to_string(graphSimplifications._bulgeLen_kAdd)));
  simplificationsParser->push_back (new OptionOneParam ("-bulge-altpath-kadd", "explore up to k+X nodes to find alternative path",  false, to_string(graphSimplifications._bulgeAltPath_kAdd))); // TODO k should not appear in that equation
  simplificationsParser->push_back (new OptionOneParam ("-bulge-altpath-covmult", "bulges of coverage <= X*cov_altpath will be removed",  false, to_string(graphSimplifications._bulgeAltPath_covMult))); 
  simplificationsParser->push_back (new OptionOneParam ("-ec-len-kmult",       "EC shorter than k*X bp are candidates to be removed",  false, to_string(graphSimplifications._ecLen_kMult)));
  simplificationsParser->push_back (new OptionOneParam ("-ec-rctc-cutoff",     "EC relative coverage coefficient (similar in spirit as tip)",  false, to_string(graphSimplifications._ecRCTCcutoff)));
  getParser()->push_back (simplificationsParser);

  IOptionsParser* graphParser = Graph::getOptionsParser(false);
  if (IOptionsParser* p = graphParser->getParser(STR_URI_INPUT))  {  p->setVisible(false); }
  if (Option* p = dynamic_cast<Option*> (graphParser->getParser(STR_KMER_ABUNDANCE_MIN)))  {  p->setDefaultValue ("2"); }
  char *hidden_env_variable_out_tmp = std::getenv("ITRAS_OUT_TMP");
  if (hidden_env_variable_out_tmp){
    if (Option* p = dynamic_cast<Option*> (graphParser->getParser(STR_URI_OUTPUT_TMP)))  {  p->setDefaultValue ((string)hidden_env_variable_out_tmp); }
  }
  getParser()->push_back(graphParser, 1);
}

void ItrAs::ParseOptions(){
  mink = getInput()->getInt(STR_MINK);
  maxk = getInput()->getInt(STR_MAXK);
  step = getInput()->getInt(STR_STEP);
  min_contig = getInput()->getInt(STR_MINCONTIG);
  read_file = getInput()->getStr(STR_URI_INPUT);
  directory = getInput()->getStr("-out");
  bool_precorrect = getInput()->getInt(STR_BOOL_PRECORRECT);
  celltype = getInput()->getStr(STR_CELLTYPE);
  matchtype = getInput()->getStr(STR_MATCHTYPE);
  aggressive = getInput()->getStr(STR_AGGRESIVENESS);
  precor_stages = getInput()->getStr(STR_PRECOR_STAGES);
  precor_trim = getInput()->getStr(STR_PRECOR_TRIM);
  MakeDir(directory);
  SetKVals();
}

struct Parameter{
  Parameter (ItrAs& itr) : itras_(itr){}
  ItrAs&         itras_;
};

template<size_t span> 
struct ItrAsFunctor  {void operator ()  (Parameter parameter){
  ItrAs& itras_ = parameter.itras_;
  itras_.getInput()->setStr("-out-dir", itras_.directory);
  itras_.getInput()->setInt("-verbose", 0);
  itras_.getInput()->setStr("-out-tmp", itras_.directory);
  itras_.getInput()->setStr(STR_TRAVERSAL_KIND, "contig");//assemble will now simplify the graph
  typedef GraphUnitigsTemplate<span> GraphType;
  GraphType graph;
  
  TIME_INFO (itras_.getTimeInfo(), "graph construction");
  int kSize = itras_.mink;
  itras_.getInput()->setStr(STR_URI_INPUT, itras_.read_file);
  itras_.getInput()->setInt(STR_KMER_SIZE, kSize);
  itras_.getInput()->setStr("-out", itras_.directory+FormatString("/graph-%d", kSize));
  graph = GraphType::create(itras_.getInput());
  string contigFile = itras_.assemble<GraphType, NodeGU, EdgeGU, span>(graph);
  uint nb_threads = 1;  // doesn't matter because for now link_tigs is single-threaded
  bool verbose = false;
  link_tigs<span>(contigFile, kSize, nb_threads, itras_.nbContigs, verbose);
  RemoveTmpFiles(itras_.directory, "");
  ///We have our original contigs now
  for(unsigned int i = 1; i < itras_.kVals.size(); ++i){
    int kSize = itras_.kVals[i];
    ///create the graph
    TIME_INFO (itras_.getTimeInfo(), "graph construction");
    itras_.getInput()->setInt(STR_KMER_SIZE, kSize);
    itras_.getInput()->setInt("-abundance-min", 0);
    ///put reads and contigs together:
    std::ifstream ifile(itras_.read_file.c_str());
    std::ofstream ofile(contigFile.c_str(), std::ios::app);
    ofile <<ifile.rdbuf();
    ifile.close(); ofile.close();
    ///
    itras_.getInput()->setStr(STR_URI_INPUT, itras_.directory+FormatString("/graph-%d.contigs.fa", itras_.kVals[i-1]));
    itras_.getInput()->setStr("-out", itras_.directory+FormatString("/graph-%d", kSize));
    graph = GraphType::create(itras_.getInput());
    ///now we need to extend our contigs
    std::string fileToDel = contigFile;
    contigFile = itras_.assemble<GraphType, NodeGU, EdgeGU, span>(graph);
    // link contigs
    uint nb_threads = 1;  // doesn't matter because for now link_tigs is single-threaded
    bool verbose = false;
    link_tigs<span>(contigFile, kSize, nb_threads, itras_.nbContigs, verbose);
    RemoveTmpFiles(itras_.directory, fileToDel);
  }
  FinalizeContigFile(itras_.directory, contigFile);
}};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void ItrAs::execute ()
{
  ParseOptions();

  if(bool_precorrect)//read precorrection
    ec_karect(read_file, celltype, matchtype, aggressive, directory, precor_stages, precor_trim);

  Integer::apply<ItrAsFunctor, Parameter> (maxk, Parameter(*this));
}


template <typename Graph_type, typename Node, typename Edge, size_t span>
void ItrAs::assembleFrom(Node startingNode, Graph_type& graph, IBank *outputBank)
{
    unsigned int isolatedCutoff = std::max(2*(unsigned int)graph.getKmerSize(), (unsigned int)150);
    bool isolatedLeft, isolatedRight;
    float coverage = 0;
    string sequence = graph.simplePathBothDirections(startingNode, isolatedLeft, isolatedRight, true, coverage);

    Sequence seq (Data::ASCII);
    seq.getData().setRef ((char*)sequence.c_str(), sequence.size());
    stringstream ss1;
    ss1 << nbContigs  << " LN:i:" << sequence.size() << " KC:i:" << (unsigned int)(coverage*(sequence.size()-graph.getKmerSize()+1)) << " km:f:" << fixed << std::setprecision(3) << coverage ;
    seq._comment = ss1.str();
    
    unsigned int lenTotal = sequence.size();
    if (lenTotal > min_contig && (lenTotal > isolatedCutoff || (lenTotal <= isolatedCutoff && (!(isolatedLeft && isolatedRight))) || keepIsolatedTigs))
    {
        outputBank->insert (seq);
        nbContigs += 1;
        totalNt   += lenTotal;
        if (lenTotal > maxContigLen)      { maxContigLen      = lenTotal;   }
    }
    else
        nbSmallContigs++;
    return;
}

template <typename Graph_type, typename Node, typename Edge, size_t span>
string ItrAs::assemble (/*const, removed because Simplifications isn't const anymore*/ Graph_type& graph)
{
    TIME_INFO (getTimeInfo(), "assembly");
    string output = (getInput()->get(STR_URI_OUTPUT) ?
        getInput()->getStr(STR_URI_OUTPUT) :
        System::file().getBaseName (getInput()->getStr(STR_URI_INPUT)) 
                    )+ ".contigs.fa";

    /** We create the output bank. Note that we could make this a little bit prettier
     *  => possibility to save the contigs in specific output format (other than fasta).  */
    ///TODO:
    IBank* outputBank = new BankFasta (output);
    LOCAL (outputBank);

    /** We set the fasta line size. */
    BankFasta::setDataLineSize (getInput()->getInt (STR_FASTA_LINE_SIZE));

    bool simplifyGraph = false;

    string traversalKind = "unitig"; // we output unitigs of the simplified graph or the original graph
    simplifyGraph = getInput()->getStr(STR_TRAVERSAL_KIND).compare("contig") == 0;

    nbContigs         = 0;
    nbSmallContigs    = 0;
    totalNt           = 0;
    maxContigLen      = 0;

    keepIsolatedTigs = getParser()->saw(STR_KEEP_ISOLATED);
    string str_tipRemoval = "", str_bubbleRemoval = "", str_ECRemoval = "";

    /** We get an iterator over all nodes . */
    ProgressGraphIteratorTemplate<Node,ProgressTimerAndSystem> itNode (graph.Graph_type::iterator(), progressFormat0);

    // if we want unitigs, then don't simplify the graph; else do it
    if (simplifyGraph){
      int nbCores = getInput()->getInt(STR_NB_CORES);
      bool verbose=false;
      Simplifications<Graph_type,Node,Edge> graphSimplifications(&graph, nbCores, verbose);
      if (getParser()->saw("-no-tip-removal"))
          graphSimplifications._doTipRemoval = false;
      if (getParser()->saw("-no-bulge-removal"))
          graphSimplifications._doBulgeRemoval = false;
      if (getParser()->saw("-no-ec-removal"))
          graphSimplifications._doECRemoval = false;
      if (getParser()->saw("-tip-len-topo-kmult"))
        graphSimplifications._tipLen_Topo_kMult = getInput()->getDouble("-tip-len-topo-kmult");
      if (getParser()->saw("-tip-len-rctc-kmult"))
        graphSimplifications._tipLen_RCTC_kMult = getInput()->getDouble("-tip-len-rctc-kmult");
      if (getParser()->saw("-tip-rctc-cutoff"))
        graphSimplifications._tipRCTCcutoff = getInput()->getDouble("-tip-rctc-cutoff");
      if (getParser()->saw("-bulge-len-kmult"))
        graphSimplifications._bulgeLen_kMult = getInput()->getDouble("-bulge-len-kmult");
      if (getParser()->saw("-bulge-len-kadd"))
        graphSimplifications._bulgeLen_kAdd = getInput()->getDouble("-bulge-len-kadd");
      if (getParser()->saw("-bulge-altpath-kadd"))
        graphSimplifications._bulgeAltPath_kAdd = getInput()->getDouble("-bulge-altpath-kadd");
      if (getParser()->saw("-bulge-altpath-covMult"))
        graphSimplifications._bulgeAltPath_covMult = getInput()->getDouble("-bulge-altpath-covMult");
      if (getParser()->saw("-ec-len-kmult"))
        graphSimplifications._ecLen_kMult = getInput()->getDouble("-ec-len-kmult");
      if (getParser()->saw("-ec-rctc-cutoff"))
        graphSimplifications._ecRCTCcutoff = getInput()->getDouble("-ec-rctc-cutoff");
                                                                                                                                                                                                                
      graphSimplifications.simplify();

      str_tipRemoval = graphSimplifications.tipRemoval;
      str_bubbleRemoval = graphSimplifications.bubbleRemoval;
      str_ECRemoval = graphSimplifications.ECRemoval;
    }

    /** We loop over all nodes. */
    for (itNode.first(); !itNode.isDone(); itNode.next()){
        Node node = itNode.item();
        if (graph.unitigIsMarked(node))  {  continue;   }
        if (graph.isNodeDeleted(node)) { continue; }
        assembleFrom<Graph_type, Node, Edge, span>(node, graph, outputBank);
    }
    
    /** We add the input parameters to the global properties. */
    // getInfo()->add (1, getInput());

    // /** We gather some statistics. */
    getInfo()->add (1, "stats");
    getInfo()->add (2, "traversal",         "%s", getInput()->getStr(STR_TRAVERSAL_KIND).c_str());
    getInfo()->add (2, "nb_contigs",         "%d", nbContigs);
    getInfo()->add (2, "nb_small_contigs_discarded","%d", nbSmallContigs);
    getInfo()->add (2, "nt_assembled",      "%ld", totalNt);
    getInfo()->add (2, "max_length",        "%d", maxContigLen);
    getInfo()->add (2, "graph simpification stats");
    getInfo()->add (3, "tips removed",          "%s", str_tipRemoval.c_str());
    getInfo()->add (3, "bulges removed",          "%s", str_bubbleRemoval.c_str());
    getInfo()->add (3, "EC removed",          "%s", str_ECRemoval.c_str());
    getInfo()->add (2, "assembly traversal stats");

    return output;
}