//arguments: reference_pdb  template_site target_pdb reference_atom
//OR
//arguments: reference_pdb template_site mutations target_pdb reference_atom

#include "newGASS-base/newGASS.hpp"
#include "newGASS-base/site.hpp"
#include <iostream>
#include <string>
#include <map>
#include <unistd.h>
#include <string.h>
#include <regex>
#include <sstream>
#include <fstream>
#include <algorithm> 
#include <cctype>
#include <locale>
#include <filesystem>
#include <boost/program_options.hpp>
#include <curl/curl.h>
#include "json.hpp"
#include <boost/thread/thread.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio.hpp>
#include <functional>

using namespace boost::program_options;
using namespace std;
namespace fs = std::filesystem;


// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
static inline std::string trim(std::string &s) {
    rtrim(s);
    ltrim(s);
    return s;
}

size_t write_data(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    size_t written;
    written = fwrite(ptr, size, nmemb, stream);
    return written;
}


bool isPdbCode(const string& code) {
    const regex pdbCodeRegex("[0-9][a-zA-Z0-9]{3}");
    return regex_match(code, pdbCodeRegex);
}

bool isUniprotCode(const string& code) {
    const regex uniprotCodeRegex("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}");
    return regex_match(code, uniprotCodeRegex);
}

string getPdbUrlFromRcsb(const string& code) {
    return "https://files.rcsb.org/download/" + code + ".pdb";
}
using json = nlohmann::json;

static size_t WriteCallback(void *contents, size_t size, size_t nmemb, void *userp) {
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
}

std::string getPdbUrlFromAlphafold(const std::string& code) {

    std::string url = "";
    std::string pdbUrl = "";

    // Build the URL for the Alphafold API
    url = "https://www.alphafold.ebi.ac.uk/api/prediction/" + code;
    // Initialize the curl handle
    CURL *curl = curl_easy_init();
    if (curl) {
        // Set the URL for the curl handle
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        // Set the callback function for the curl handle to capture the response
        std::string response_string;
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response_string);
        // Perform the curl request
        CURLcode res = curl_easy_perform(curl);
        if (res != CURLE_OK) {
            std::cerr << "curl_easy_perform() failed: " << curl_easy_strerror(res) << std::endl;
        }
        else {
            // Parse the JSON response
            json response_json = json::parse(response_string);
            if (response_json.size() > 0) {
                pdbUrl = response_json[0]["pdbUrl"].get<std::string>();
            }
        }
        // Clean up the curl handle
        curl_easy_cleanup(curl);
    }
    else {
        std::cerr << "Could not initialize curl" << std::endl;
    }

    return pdbUrl;
}
void downloadPdb(const string& code, const string& pathToSave) {
    string url = "";
    if (isPdbCode(code)) {
        url = getPdbUrlFromRcsb(code);
    } else if (isUniprotCode(code)) {
        url = getPdbUrlFromAlphafold(code);
    } else {
        cout << "Format error" << endl;
        exit(1);
    }

    CURL *curl = curl_easy_init();
    if (curl) {
        FILE *fp = fopen(pathToSave.c_str(), "wb");
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
        CURLcode res = curl_easy_perform(curl);
        fclose(fp);
        if (res != CURLE_OK) {
            cout << "Download error: " << curl_easy_strerror(res) << endl;
            exit(1);
        }
        curl_easy_cleanup(curl);
    } else {
        cout << "Unable to initialize CURL" << endl;
        exit(1);
    }
}
void transfer(char *linha,AtomoCompat &aux){
    aux.cadeia=linha[21]; //Adicionando a cadeia
    //Pegando o nome do atomo;
    char temp5[4];
    int j=0;
    for(int i=12;i<16;i++){
        if(linha[i]!=' '){
            temp5[j]=linha[i];
            j++;
        }
    }
    temp5[j]='\0';
    strcpy(aux.atomo,temp5);

    // Pegando o nome do aminoacido;
    char temp[4];
    j=0;
    for(int i=17;i<20;i++){
        if(linha[i]!=' '){
            temp[j]=linha[i];
            j++;
        }
    }
    temp[j]='\0';
    strcpy(aux.amino,temp);

    // Coordenada X;
    char temp2[8];
    j=0;
    for(int i=30;i<38;i++){
        if(linha[i]!=' '){
            temp2[j]=linha[i];
            j++;
        }
    }
    temp2[j]='\0';
    float X = atof(temp2);
    aux.x=X;

    //Coordenada Y
    char temp3[8];
    j=0;
    for(int i=38;i<46;i++){
        if(linha[i]!=' '){
            temp3[j]=linha[i];
            j++;
        }
    }
    temp3[j]='\0';
    float Y = atof(temp3);
    aux.y=Y;

    // Coordenada Z
    char temp4[8];
    j=0;
    for(int i=46;i<54;i++){
        if(linha[i]!=' '){
            temp4[j]=linha[i];
            j++;
        }
    }
    temp4[j]='\0';
    float Z = atof(temp4);
    aux.z=Z;

    // ID do átomo
    char temp6[4];
    j=0;
    for(int i=22;i<26;i++){
        if(linha[i]!=' '){
            temp6[j]=linha[i];
            j++;
        }
    }
    temp6[j]='\0';
    int ID = atoi(temp6);
    aux.atomo_ID=ID;
};


void generateDat(string txtFilePath, string datFilePath){
    char nome_proteina[10], nome_arquivo[10];
	char nometemp [6];

    char buff[10];  // variável para receber cada linha do arquivo
    int nlinha = 0; // conta o numero de linhas
    AtomoCompat aux;

	
    ifstream inFileProteina;
    inFileProteina.open(txtFilePath);
    if(!inFileProteina){
        cout<<"Erro ao abrir o arquivo: "<<txtFilePath<<endl;
        exit(1);
    }
    else{
        char linha[100];
        ofstream outFile(datFilePath,ios::binary);
        int countline = 0;
        while(!inFileProteina.eof()){
            inFileProteina.getline(linha,100);
            if(countline<=3){
                outFile.write(reinterpret_cast<const char *>(&linha), sizeof(linha));
                countline++;

                } else{
                    transfer(linha,aux);
                    countline++;
                    outFile.write(reinterpret_cast<const char *>(&aux), sizeof(AtomoCompat));
                }
        }
        outFile.write(reinterpret_cast<const char *>(&countline), sizeof(int));
        outFile.close();
        inFileProteina.close();
    }
}

void generateTxt(string pdbFilePath,string txtFilePath, string referenceAtom){
    const map<string,string> dict_lha = {{"ALA","CB"},{"ARG","CZ"},{"ASN","CG"},
                                        {"ASP","CG"}, {"CYS","SG"}, {"GLN","CD"},
                                        {"GLU","CD"}, {"GLY","CA"}, {"HIS","CE1"},
                                        {"ILE","CD1"},{"LEU","CD1"}, {"LYS","NZ"},
                                        {"MET","CE"}, {"PHE","CZ"}, {"PRO","CG"},
                                        {"SER","OG"},{"THR","CG2"}, {"TRP","CH2"},
                                        {"TYR","OH"}, {"VAL","CG1"}};
    ifstream arqPdb(pdbFilePath);
    ofstream arqTxt(txtFilePath);

    string proteinName = "NULL";
    string ecNumber = "NULL";
    string unp = "NULL";
    string resolution = "NULL";

    bool first = true;
    while(!arqPdb.eof()){
        string line;
        getline(arqPdb,line);
        if(line.find("HEADER")!= string::npos)
            proteinName = line.substr(62,4);
        if(line.find("EC:")!= string::npos && line.find("COMPND")!= string::npos){
            ecNumber = line.substr(15);
            trim(ecNumber);
            while(ecNumber.find(";")!= string::npos)
                ecNumber.replace(ecNumber.find(";"),1," ");
        }
        if(line.find("UNP")!= string::npos && line.find("DBREF")!= string::npos)
            unp = line.substr(33,6);
        if(line.find("RESOLUTION.")!= string::npos && line.find("REMARK")!= string::npos){
            resolution = line.substr(26,5);
            trim(resolution);
        }

        if(line.substr(0,4) == "ATOM" ){
            if(first){
                arqTxt << proteinName << "\n" << ecNumber << "\n" << unp << "\n" << resolution << endl;
                first = false;
            }
            string res = line.substr(17,3);
            string ra;
            if(referenceAtom=="CA")
                ra = "CA";
            else
                ra = dict_lha.at(res);
            string atomTipe = line.substr(12,4);
            trim(atomTipe);
            if(atomTipe== ra)
                arqTxt << line << endl;
        }
        

        
    }
}

site readTemplate(string templateSite, string referencePdbFileName){
    site temp;
    string word;
    stringstream ss(templateSite);
    int i=0;
    Atomo aux;

    while(ss.good()){
        if(i%3==0){
            getline(ss,word,',');
            aux.amino=word;
        }
        else if(i%3==1){
            getline(ss,word,',');
            aux.atomo_ID=stoi(word);
        }
        else if(i%3==2){
            getline(ss,word,';');
            aux.cadeia=word[0];
            temp.residuos.push_back(aux);
        }
        i++;
    }
    temp.size=temp.residuos.size();
    GASS::readTemplateFile(referencePdbFileName,temp);
    return temp;
}

vector<site> readTemplates(string jsonFilePath, string cacheFolder) {
    vector<site> templates;
    ifstream jsonFile(jsonFilePath);
    if (!jsonFile) {
        cerr << "Failed to open JSON file: " << jsonFilePath << endl;
        exit(1);
    }
    json jsonData;
    try {
        jsonFile >> jsonData;
    }
    catch (const std::exception& e) {
        cerr << "Failed to parse JSON file: " << e.what() << endl;
        exit(1);
    }
    if (jsonData.find("TEMPLATES") != jsonData.end()) {
        for (auto& templateJson : jsonData["TEMPLATES"]) {
            site templateObj;
            if (templateJson.find("REFERENCE_PDB") != templateJson.end()) {
                templateObj.pdbId = templateJson["REFERENCE_PDB"].get<std::string>();
            }
            if (templateJson.find("RESIDUES") != templateJson.end()) {
                for (auto& residueJson : templateJson["RESIDUES"]) {
                    Atomo residueObj;
                    if (residueJson.find("AMINO") != residueJson.end()) {
                        residueObj.amino = residueJson["AMINO"].get<std::string>();
                    }
                    if (residueJson.find("ID") != residueJson.end()) {
                        residueObj.atomo_ID = residueJson["ID"].get<int>();
                    }
                    if (residueJson.find("CHAIN") != residueJson.end()) {
                        residueObj.cadeia = residueJson["CHAIN"].get<string>()[0];
                    }
                    templateObj.residuos.push_back(residueObj);
                }
            }
            templateObj.size = templateObj.residuos.size();
            string referencePdbFileName = cacheFolder + "/" + templateObj.pdbId + "/protein.dat";
            if(!fs::exists(cacheFolder+templateObj.pdbId)){
                fs::create_directory(cacheFolder+templateObj.pdbId);
                downloadPdb(templateObj.pdbId, cacheFolder+templateObj.pdbId+"/protein.pdb");
                generateTxt(cacheFolder+templateObj.pdbId+"/protein.pdb",cacheFolder+templateObj.pdbId+"/protein.txt","CA");
                generateDat(cacheFolder+templateObj.pdbId+"/protein.txt",cacheFolder+templateObj.pdbId+"/protein.dat");
            }
            GASS::readTemplateFile(referencePdbFileName,templateObj);
            templates.push_back(templateObj);
        }
    }
    return templates;
}
string getValueFromConfig(string configFilePath, string key) {
    ifstream inFile(configFilePath);
    if(!inFile){
        cout<<"Erro ao abrir o arquivo"<<configFilePath<<endl;
        exit(1);
    }
    else{
        json j;
        inFile >> j;

       return j[key].get<string>();
    }
    return "";
}

GASS::Parameters readConfigParameters(string configFilePath){
    GASS::Parameters param;

    ifstream inFile(configFilePath);
    if(!inFile){
        cout<<"Erro ao abrir o arquivo: "<<configFilePath<<endl;
        exit(1);
    }
    else{
        json j;
        inFile >> j;

        param.AG_NUMBER_OF_GENERATIONS = j["AG"]["NUMBER_OF_GENERATIONS"].get<int>();
        param.AG_POPULATION_SIZE = j["AG"]["POPULATION_SIZE"].get<int>();
        param.AG_MUTATION_RATE = j["AG"]["MUTATION_RATE"].get<float>();
        param.AG_CROSSOVER_RATE = j["AG"]["CROSSOVER_RATE"].get<float>();
        param.AG_NUMBER_OF_ELITE_INDIVIDUALS = j["AG"]["NUMBER_OF_ELITE_INDIVIDUALS"].get<int>();
        param.AG_TOURNAMENT_SIZE = j["AG"]["TOURNAMENT_SIZE"].get<int>();
    }

    return param;
}

GASS::Parameters getDefaultParameters(){
    return readConfigParameters("config.json");
}

GASS::Parameters getDefaultParameters(string configFilePath){
    return readConfigParameters(configFilePath);
}

void run(site temp, Repositorio rep, string outputFilePath){
    std::set<Individuo> results;
    ofstream ofs(outputFilePath);

    GASS::runOneToOne(&temp, &rep, &results, getDefaultParameters());

    int contador=1;
    int tam_sitiof = temp.residuos.size();
    int i=0;
    for (Individuo ind:results)
    {
        ofs << contador << "\t" << i << "\t";
        ofs << std::setprecision(3) << ind.getFitness() << "\t";
        if (ind.getFitness() != 1000)
        {
            ofs << ind.cromossomo[0].amino << " " << ind.cromossomo[0].atomo_ID << " " << ind.cromossomo[0].cadeia;
            for (int j = 1; j < tam_sitiof; j++)
                ofs << ";" << ind.cromossomo[j].amino << " " << ind.cromossomo[j].atomo_ID << " " << ind.cromossomo[j].cadeia;
            ofs << "\t" << temp.pdbId << "\t";
            ofs << temp.residuos[0].amino << " " << temp.residuos[0].atomo_ID << " " << temp.residuos[0].cadeia;
            for (int j = 1; j < tam_sitiof; j++)
                ofs << ";" << temp.residuos[j].amino << " " << temp.residuos[j].atomo_ID << " " << temp.residuos[j].cadeia;
        }
        else
        {
            ofs << " - "
                << " "
                << " - "
                << " "
                << " - ";
            for (int j = 1; j < tam_sitiof; j++)
                ofs << ";"
                    << " - "
                    << " "
                    << " - "
                    << " "
                    << " - ";
            ofs << "\t" << temp.pdbId << "\t";
            ofs << temp.residuos[0].amino << " " << temp.residuos[0].atomo_ID << " " << temp.residuos[0].cadeia;
            for (int j = 1; j < tam_sitiof; j++)
                ofs << ";" << temp.residuos[j].amino << " " << temp.residuos[j].atomo_ID << " " << temp.residuos[j].cadeia;
        
        }
        //ofs << "\t" << ecnumber << "\t" << uniprot << "\t" << resolution << "\t" << tam_sitiof;
        ofs << std::endl;
        i++;
    }
}

void run(vector<site> temps, Repositorio rep, string outputFilePath, GASS::Parameters param, int nThreads){
    vector<std::set<Individuo>> results (temps.size());
    ofstream ofs(outputFilePath);
    boost::asio::thread_pool pool(nThreads);
    auto task = [&] (int i) {
        GASS::runOneToOne(&temps[i], &rep, &results[i], param);
    };
    for(int i=0;i<temps.size();i++){
        boost::asio::post(pool, std::bind(task, i));
    }
    pool.join();
        
    int contador=1;
        
    for(int j=0;j<temps.size();j++){
        site temp = temps[j];
        int tam_sitiof = temp.residuos.size();
        int i=0;
        for (Individuo ind:results[j])
        {
            ofs << contador << "\t" << i << "\t";
            ofs << std::setprecision(3) << ind.getFitness() << "\t";
            if (ind.getFitness() != 1000)
            {
                ofs << ind.cromossomo[0].amino << " " << ind.cromossomo[0].atomo_ID << " " << ind.cromossomo[0].cadeia;
                for (int j = 1; j < tam_sitiof; j++)
                    ofs << ";" << ind.cromossomo[j].amino << " " << ind.cromossomo[j].atomo_ID << " " << ind.cromossomo[j].cadeia;
                ofs << "\t" << temp.pdbId << "\t";
                ofs << temp.residuos[0].amino << " " << temp.residuos[0].atomo_ID << " " << temp.residuos[0].cadeia;
                for (int j = 1; j < tam_sitiof; j++)
                    ofs << ";" << temp.residuos[j].amino << " " << temp.residuos[j].atomo_ID << " " << temp.residuos[j].cadeia;
            }
            else
            {
                ofs << " - "
                    << " "
                    << " - "
                    << " "
                    << " - ";
                for (int j = 1; j < tam_sitiof; j++)
                    ofs << ";"
                        << " - "
                        << " "
                        << " - "
                        << " "
                        << " - ";
                ofs << "\t" << temp.pdbId << "\t";
                ofs << temp.residuos[0].amino << " " << temp.residuos[0].atomo_ID << " " << temp.residuos[0].cadeia;
                for (int j = 1; j < tam_sitiof; j++)
                    ofs << ";" << temp.residuos[j].amino << " " << temp.residuos[j].atomo_ID << " " << temp.residuos[j].cadeia;
            
            }
            //ofs << "\t" << ecnumber << "\t" << uniprot << "\t" << resolution << "\t" << tam_sitiof;
            ofs << std::endl;
            i++;
        }
        contador++;
    }
    results.clear();
    
    ofs.close();
}


void readMutations(string mutations, site& temp){
    const map<string, int> AMINO_CODES_MAP = {
            {"ALA", 0}, {"ARG", 1}, {"ASN", 2}, {"ASP", 3}, {"CYS", 4}, {"GLN", 5}, {"GLU", 6}, {"GLY", 7}, {"HIS", 8}, {"ILE", 9}, {"LEU", 10}, {"LYS", 11}, {"MET", 12}, {"PHE", 13}, {"PRO", 14}, {"SER", 15}, {"THR", 16}, {"TRP", 17}, {"TYR", 18}, {"VAL", 19}
        };
    stringstream ss(mutations);
    while(ss.good()){
        string res1;
        string res2;
        getline(ss, res1, ',');
        if(AMINO_CODES_MAP.find(res1.substr(0,3)) != AMINO_CODES_MAP.end()){
            getline(ss, res2, ';');
            if(AMINO_CODES_MAP.find(res1.substr(0,3)) != AMINO_CODES_MAP.end())
                temp.substitutions[AMINO_CODES_MAP.at(res1.substr(0,3))].push_back(AMINO_CODES_MAP.at(res2.substr(0,3)));
        }
    }
}

void printTemplateInfo(site temp){
 //show template infos
    cout << "Template: " << temp.pdbId << endl;
    cout << "Residues: " << temp.residuos.size() << endl;
    for(int i=0; i<temp.residuos.size(); i++){
        cout << temp.residuos[i].amino << " " << temp.residuos[i].atomo_ID << " " << temp.residuos[i].cadeia << endl;
    }
    cout << "Substitutions: " << temp.substitutions.size() << endl;
    for(int i=0; i<temp.substitutions.size(); i++){
        cout << i << " -> ";
        for(int j=0; j<temp.substitutions[i].size(); j++){
            cout << temp.substitutions[i][j] << " ";
        }
        cout << endl;
    }
}

void printTargetInfo(Repositorio rep){
    const map<string, int> AMINO_CODES_MAP = {
            {"ALA", 0}, {"ARG", 1}, {"ASN", 2}, {"ASP", 3}, {"CYS", 4}, {"GLN", 5}, {"GLU", 6}, {"GLY", 7}, {"HIS", 8}, {"ILE", 9}, {"LEU", 10}, {"LYS", 11}, {"MET", 12}, {"PHE", 13}, {"PRO", 14}, {"SER", 15}, {"THR", 16}, {"TRP", 17}, {"TYR", 18}, {"VAL", 19}
        };
    //show target infos
    cout << "Target: " << rep.pdbId << endl;
    cout << "Residues: " << endl;
    for(int i=0; i<rep.rep[0].size(); i++){
        cout << i << ": " << rep.rep[0][i].size() << endl;
    }
}


int main(int argc, char* argv[]){
    std::string appFolder ="";
    if(!fs::exists(fs::current_path()/"config.json")){
        std::cout<<"config file not found on path: "<<fs::current_path()<<std::endl;
        std::cout<<"using default path: ~/.gasscli/config"<<std::endl;
        std::string username = getlogin();
        appFolder = "/home/" + username + "/.gasscli/";
    }
    else{
        appFolder = fs::current_path().string() + "/";
    }
    
    if(!fs::exists(appFolder)){
        throw std::runtime_error("App folder not found on path: " + appFolder);
    }
    options_description run_opt("run options");
    run_opt.add_options()
        ("help,h", "produce help message")
        ("cache,c", value<string>(), "set cache folder, default is .gasscli/cache/")
        ("reference_pdb,r", value<string>(), "set reference pdb")
        ("template_site,s", value<string>(), "set template site")
        ("mutations", value<string>(), "set mutations")
        ("target_pdb,t", value<string>(), "set target pdb")
        ("reference_atom", value<string>(), "set reference atom")
        ("output,o", value<string>(), "set output file, default is output.txt")
        ("config_file", value<string>(), "set the path to the config json file to be used")
        ("num_threads,j", value<int>(), "set the number of threads to be used, default is 1")
    ;
    options_description download_opt("download options");
    download_opt.add_options()
        ("help,h", "produce help message")
        ("cache_folder,c", value<string>(), "set cache folder to download to, default is .gasscli/cache/")
        ("pdb_list,p", value< vector<string> >()->multitoken(), "set pdb list")
        ("reference_atom", value<string>(), "set reference atom, default is CA")
        ("file_list,f", value<string>(), "set the path to a file containing a list of pdb files to download")
        ("force", "force download even if .pdb file already exists")

    ;
    options_description prepare_opt("prepare options");
    prepare_opt.add_options()
        ("help,h", "produce help message")
        ("cache_folder,c", value<string>(), "set cache folder to add the .dat files to, default is .gasscli/cache/")
        ("pdb_code,p", value<string>(), "set pdb code")
        ("reference_atom", value<string>(), "set reference atom, default is CA")
        ("pdb_file_path,f", value<string>(), "set the path to the .pdb file to be prepared")
        ("force", "force prepare even if .dat file already exists")
    ;

    if(argc>0 && string(argv[1])=="run"){
        
        string configFilePath = "";
        string referencePdb = "";
        string templateSite = "";
        string mutations = "";
        string targetPdb = "";
        string referenceAtom = "CA";
        string cacheFolder = "";
        string outputFilePath = "";
        int numThreads = 1;

        variables_map vm;
        store(parse_command_line(argc, argv, run_opt), vm);
        notify(vm);

        if (vm.count("help")) {
            cout << run_opt << "\n";
            return 1;
        }
        if(vm.count("num_threads")){
            numThreads = vm["num_threads"].as<int>();
        }
        if(vm.count("config_file")){
            configFilePath = vm["config_file"].as<string>();
        }
        if(vm.count("cache")){
            cacheFolder = vm["cache"].as<string>();
            if(cacheFolder[cacheFolder.size()-1] != '/'){
                cacheFolder += "/";
            }
        }
        if(vm.count("reference_pdb")){
            referencePdb = vm["reference_pdb"].as<string>();
        }
        if(vm.count("template_site")){
            templateSite = vm["template_site"].as<string>();
        }
        if(vm.count("mutations")){
            mutations = vm["mutations"].as<string>();
        }
        if(vm.count("target_pdb")){
            targetPdb = vm["target_pdb"].as<string>();
        }
        if(vm.count("reference_atom")){
            referenceAtom = vm["reference_atom"].as<string>();
        }
        if(vm.count("output")){
            outputFilePath = vm["output"].as<string>();
        }

        if(configFilePath == "" && (referencePdb == "" || templateSite == "" || targetPdb == "")){
            cout << "Missing arguments" << endl;
            cout << run_opt << "\n";
            return 1;
        }

        if(configFilePath == ""){
            cout << "No config file specified" << endl;
            cout << "Using default config file" << endl;
            configFilePath = appFolder + "config.json";
            cout<<configFilePath<<"\n";
        }
        if(cacheFolder == ""){
            cacheFolder = getValueFromConfig(configFilePath, "DEFAULT_CACHE_FOLDER");
        }
        if(!fs::exists(cacheFolder)){
            fs::create_directory(cacheFolder);
        }
        vector<site> templates;
        if(referencePdb == ""){
            templates =readTemplates(configFilePath, cacheFolder);
        }
        else{
            if(!fs::exists(cacheFolder+referencePdb)){
                fs::create_directory(cacheFolder+referencePdb);
                downloadPdb(referencePdb, cacheFolder+referencePdb+"/protein.pdb");
                generateTxt(cacheFolder+referencePdb+"/protein.pdb",cacheFolder+referencePdb+"/protein.txt",referenceAtom);
                generateDat(cacheFolder+referencePdb+"/protein.txt",cacheFolder+referencePdb+"/protein.dat");
            }
            templates.push_back(readTemplate(templateSite, cacheFolder+referencePdb+"/protein.dat"));
        }
        if (targetPdb == ""){
            targetPdb = getValueFromConfig(configFilePath, "TARGET_PDB");
        }
        if(!fs::exists(cacheFolder+targetPdb)){
            fs::create_directory(cacheFolder+targetPdb);
            downloadPdb(targetPdb, cacheFolder+targetPdb+"/protein.pdb");
            generateTxt(cacheFolder+targetPdb+"/protein.pdb",cacheFolder+targetPdb+"/protein.txt",referenceAtom);
            generateDat(cacheFolder+targetPdb+"/protein.txt",cacheFolder+targetPdb+"/protein.dat");
        }
        
        if(mutations!="")
            readMutations(mutations, templates[0]);
    
        //printTemplateInfo(temp);

        Repositorio repositorio;
        repositorio.readRepository(cacheFolder+targetPdb+"/protein.dat");

        if(outputFilePath == ""){
            if(configFilePath != "")
                outputFilePath = getValueFromConfig(configFilePath, "DEFAULT_RESULT_FILE");
            else
                outputFilePath = "output.txt";
        }
        //printTargetInfo(repositorio);
        if(templates.size() == 1)
            run(templates[0],repositorio,outputFilePath);
        else
            ::run(templates,repositorio,outputFilePath, getDefaultParameters(configFilePath), numThreads);
    }
    else if(argc>0 && string(argv[1])=="download"){
        string cacheFolder = appFolder + "cache/";
        vector<string> pdbList;
        string referenceAtom = "CA";
        string file_list = "";
        bool force = false;

        variables_map vm;
        store(parse_command_line(argc, argv, download_opt), vm);
        notify(vm);

        if (vm.count("help")) {
            cout << download_opt << "\n";
            return 1;
        }
        if(vm.count("cache_folder")){
            cacheFolder = vm["cache_folder"].as<string>();
            if(cacheFolder[cacheFolder.size()-1] != '/'){
                cacheFolder += "/";
            }
        }
        if(vm.count("pdb_list")){
            pdbList = vm["pdb_list"].as<vector<string>>();
        }
        if(vm.count("reference_atom")){
            referenceAtom = vm["reference_atom"].as<string>();
        }
        if(vm.count("file_list")){
            file_list = vm["file_list"].as<string>();
        }
        if(vm.count("force")){
            force = true;
        }

        if(pdbList.empty() && file_list == ""){
            cout << "Missing arguments" << endl;
            cout << download_opt << "\n";
            return 1;
        }

        if(!fs::exists(cacheFolder)){
            fs::create_directory(cacheFolder);
        }
        if(!pdbList.empty()){
            for(string str : pdbList)
            {
                if(!fs::exists(cacheFolder+str)){
                    fs::create_directory(cacheFolder+str);
                    downloadPdb(str, cacheFolder+str+"/protein.pdb");
                    generateTxt(cacheFolder+str+"/protein.pdb",cacheFolder+str+"/protein.txt",referenceAtom);
                    generateDat(cacheFolder+str+"/protein.txt",cacheFolder+str+"/protein.dat");
                }
                else if(force){
                    downloadPdb(str, cacheFolder+str+"/protein.pdb");
                    generateTxt(cacheFolder+str+"/protein.pdb",cacheFolder+str+"/protein.txt",referenceAtom);
                    generateDat(cacheFolder+str+"/protein.txt",cacheFolder+str+"/protein.dat");
                }
            }
        }
        if(file_list != ""){
            ifstream file(file_list);
            string str; 
            while (std::getline(file, str))
            {
                if(!fs::exists(cacheFolder+str)){
                    fs::create_directory(cacheFolder+str);
                    downloadPdb(str, cacheFolder+str+"/protein.pdb");
                    generateTxt(cacheFolder+str+"/protein.pdb",cacheFolder+str+"/protein.txt",referenceAtom);
                    generateDat(cacheFolder+str+"/protein.txt",cacheFolder+str+"/protein.dat");
                }
                else if(force){
                    downloadPdb(str, cacheFolder+str+"/protein.pdb");
                    generateTxt(cacheFolder+str+"/protein.pdb",cacheFolder+str+"/protein.txt",referenceAtom);
                    generateDat(cacheFolder+str+"/protein.txt",cacheFolder+str+"/protein.dat");
                }
            }
        }
    }
    else if(argc>0 && string(argv[1])=="prepare"){
        string cacheFolder = appFolder + "cache/";
        string pdbCode = "";
        string referenceAtom = "CA";
        string pdb_file_path = "";
        bool force = false;

        variables_map vm;
        store(parse_command_line(argc, argv, prepare_opt), vm);
        notify(vm);

        if (vm.count("help")) {
            cout << prepare_opt << "\n";
            return 1;
        }
        if(vm.count("cache_folder")){
            cacheFolder = vm["cache_folder"].as<string>();
            if(cacheFolder[cacheFolder.size()-1] != '/'){
                cacheFolder += "/";
            }
        }
        if(vm.count("pdb_code")){
            pdbCode = vm["pdb_code"].as<string>();
        }
        if(vm.count("reference_atom")){
            referenceAtom = vm["reference_atom"].as<string>();
        }
        if(vm.count("pdb_file_path")){
            pdb_file_path = vm["pdb_file_path"].as<string>();
        }
        if(vm.count("force")){
            force = true;
        }

        if(pdbCode == "" && pdb_file_path == ""){
            cout << "Missing arguments" << endl;
            cout << prepare_opt << "\n";
            return 1;
        }

        if(!fs::exists(cacheFolder)){
            fs::create_directory(cacheFolder);
        }
        if(pdbCode != ""){
            if(!fs::exists(cacheFolder+pdbCode)){
                fs::create_directory(cacheFolder+pdbCode);
                generateTxt(pdb_file_path,cacheFolder+pdbCode+"/protein.txt",referenceAtom);
                generateDat(cacheFolder+pdbCode+"/protein.txt",cacheFolder+pdbCode+"/protein.dat");
            }
            else if(force){
                generateTxt(pdb_file_path,cacheFolder+pdbCode+"/protein.txt",referenceAtom);
                generateDat(cacheFolder+pdbCode+"/protein.txt",cacheFolder+pdbCode+"/protein.dat");
            }
        }

    }
    else{
        cout << "Usage: gasscli [instruction] [options]\n\n"
         << "Allowed instructions:\n"
         << run_opt << "\n"
         << download_opt << "\n"
         << prepare_opt << endl;
    }
    return 0;
}