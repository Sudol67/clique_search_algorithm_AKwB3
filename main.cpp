#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <utility>
#include <unordered_map>
#include <set>

using namespace std;

vector<string> sequences; //oryginalne sekwencje wczytane z pliku
vector<string> sequencesQual;
vector<string> trimmedSequences;
vector<string> trimmedQualities;
vector<vector<int>> trimmedPositions;
vector<pair<int, pair<int, string>>> windows;
unordered_map<string, vector<int>> graph;
int windowSize;


vector<string> readSequences(string filepath) {
    string line;
    string sequence;
    bool isSequence = false;
    ifstream file(filepath);
    if (file.is_open()) {
        while (getline(file, line)) {
            if (line[0] == '>') { //Jesli linia zaczyna się od znaku '>' to nie jest wczytywana
                isSequence = false;
            } else if (!isSequence) {
                isSequence = true;
                sequence = line;
                sequences.push_back(sequence);
            }
        }
        file.close();
    } else {
        cout << "Unable to open file";
    }
    return sequences;
}

vector<string> readQuality(string filepath) {
    string line;
    string quality;
    bool isSequence = false;
    ifstream file(filepath);
    if (file.is_open()) {
        while (getline(file, line)) {
            if (line[0] == '>') {  //Jesli linia zaczyna się od znaku '>' to nie jest wczytywana
                isSequence = false;
            } else if (!isSequence) {
                isSequence = true;
                quality = line;
                sequencesQual.push_back(quality);
            }
        }
        file.close();
    } else {
        cout << "Unable to open file";
    }
    return sequencesQual;
}

void trimQuality(vector<string> &sequences, vector<string> &sequencesQual, int minQuality, vector<vector<int>> &trimmedPositions) {
    for (int i = 0; i < sequences.size(); i++) {
        string trimmedSeq, trimmedQuality;
        stringstream qualStream(sequencesQual[i]);
        vector<int> positions;
        for (int j = 0; j < sequences[i].length(); j++) {
            int quality;
            qualStream >> quality;
            if (quality >= minQuality) {
                trimmedSeq += sequences[i][j];
                trimmedQuality += to_string(quality) + " ";
                positions.push_back(j+1);
            }
        }
        trimmedSequences.push_back(trimmedSeq);
        trimmedQualities.push_back(trimmedQuality);
        trimmedPositions.push_back(positions);
    }
}

vector<pair<int, pair<int, string>>> getWindows(vector<string> &trimmedSequences, vector<vector<int>> &trimmedPositions, int windowSize) {
    for (int i = 0; i < trimmedSequences.size(); i++) { //iteracja przez sekwencje. Pierwsza pętla iteruje przez wszystkie sekwencje a druga przez konkretną sekwencję w poszukiwaniu okna
        for (int j = 0; j < trimmedSequences[i].size() - windowSize + 1; j++) { //'j < trimmedSequences[i].size() - windowSize + 1' zapewnia odpowiednią długość okna któa byłą wcześniej zdefiniowana
            string window = trimmedSequences[i].substr(j, windowSize); //pobieranie okna sekwencji. Zaczynamy od pozycji 'j' w sekwencji trimmedSequences[i]. Funkcja substr przyjmuje dwa argumenty: 'j' czyli pozycja od któej zaczynamy i windowSize czyli ilość znaków które chcemy pobrać
            windows.push_back({i, {trimmedPositions[i][j], window}}); //zapisywanie okna. 'i' to numer sekwencji, trimmedPositions to pozycja rozpoczynająca okno i Window to string z sekwencją okna
        }
    }
    return windows;
}

void createGraph() {
    for (int i = 0; i < windows.size(); i++) {//pierwsza pętla iteruje po wszystkich oknach sekwencji (indeksach okien w wektorze windows)
        for (int j = i + 1; j < windows.size(); j++) { //druga pętla pozwala na porównanie każdego okna z pozostałymi
            if(windows[i].second.second != windows[j].second.second){
                continue;
            }
            if (windows[i].first == windows[j].first){ //sprawdzenie czy okna pochodzą z tej samej sekwencji. Jeśli pochodzą, wykonywane jest continue i przerywana jest dana iteracja
                continue;
            }
            int pos1 = windows[i].second.first; //porównanie odległości. Zawsze mniejsza wartość odejmowana jest od większej
            int pos2 = windows[j].second.first;
            if (pos1 > pos2) {
                swap(pos1, pos2); // zamiana miejscami, jeżeli pos1 jest większe od pos2
                }
            if (pos2 - pos1 > 10 * windowSize){
                continue; //Jeśli różnica jest większa niż dziesięciokrotność długości okna, wykonywane jest continue i przerywana jest dana iteracja
            }

            if(graph[windows[i].second.second].empty()){
                graph[windows[i].second.second].push_back(i);
            }
            graph[windows[i].second.second].push_back(j); //Jeśli wszystkie warunki zostaną spełnione, tworzone jest połączenie. Kluczem mapy 'graph' jest sekwencja okna a wartością jest wektor przechowujący indeksy okien z którymi jest połączenie
        }
    }
}

void find_star(unordered_map<string, vector<int>>& graph, vector<pair<int, pair<int, string>>>& windows) {
    int pos1, pos2;
    bool struktura = false;
    bool isStar = true;
    for (auto const& vertex : graph) {
        if (vertex.second.size() == 11) {
            struktura = true;
            string key = vertex.first;
            auto it = graph.find(key);
            if(it != graph.end()){ //najpierw porownanie pierwszego wierzcholka do wszystkich a pozniej po kolei kazdy do kazdego (oprocz 1)
                for(int i = 1;i < 5;i++){
                    pos1 = windows[it->second.at(0)].second.first;
                    pos2 = windows[it->second.at(i)].second.first;
                    if (pos1 > pos2) {
                        swap(pos1, pos2); // zamiana miejscami, jeżeli pos1 jest większe od pos2
                    }
                    if(windows[it->second.at(0)].first == windows[it->second.at(i)].first || pos2 - pos1 > 10 * windowSize){
                        isStar = false;
                    }
                }
                for(int i = 5;i < 8;i++){
                    pos1 = windows[it->second.at(1)].second.first;
                    pos2 = windows[it->second.at(i)].second.first;
                    if (pos1 > pos2) {
                        swap(pos1, pos2); // zamiana miejscami, jeżeli pos1 jest większe od pos2
                    }
                    if(windows[it->second.at(1)].first == windows[it->second.at(i)].first || pos2 - pos1 > 10 * windowSize){
                        isStar = false;
                    }
                }
                for(int i = 8;i < 10;i++){
                    pos1 = windows[it->second.at(5)].second.first;
                    pos2 = windows[it->second.at(i)].second.first;
                    if (pos1 > pos2) {
                        swap(pos1, pos2); // zamiana miejscami, jeżeli pos1 jest większe od pos2
                    }
                    if(windows[it->second.at(5)].first == windows[it->second.at(i)].first || pos2 - pos1 > 10 * windowSize){
                        isStar = false;
                    }
               }
                pos1 = windows[it->second.at(8)].second.first;
                pos2 = windows[it->second.at(10)].second.first;
                if (pos1 > pos2) {
                    swap(pos1, pos2); // zamiana miejscami, jeżeli pos1 jest większe od pos2
                }
                if(windows[it->second.at(8)].first == windows[it->second.at(10)].first || pos2 - pos1 > 10 * windowSize){
                    isStar = false;
                }
            }
            if(isStar == true){
                cout<<"Motyw: "<<key<<endl;
                for(int m = 0;m < 5;m++){
                    cout<<"W sekwencji: "<<windows[it->second.at(m)].first+1<<" znaleziony na pozycji: "<<windows[it->second.at(m)].second.first<<endl;
                }
                break;
            }
            else{
                cout << "Nie znaleziono struktury typu gwiazda." << endl;
            }
        }

    }
    if(struktura == false){
        cout << "Brak struktury typu gwiazda. Sprobuj uzyc innych ustawien lub wczytaj inny plik" << endl;
    }
}

int main()
{
    vector<string> sequences = readSequences("C:/Users/Basia/Desktop/AKwB_3/Gotowe/Instancja2_seq.fasta");
    vector<string> sequencesQual = readQuality("C:/Users/Basia/Desktop/AKwB_3/Gotowe/Instancja2_jakosc.qual");

    int minQuality;
    cout << "Podaj minimalny prog jakosci nukleotydu: ";
    cin >> minQuality;

    trimQuality(sequences, sequencesQual, minQuality, trimmedPositions);

    int tmp = 0;
    while(tmp == 0){
        cout << "Podaj rozmiar okna: ";
        cin >> windowSize;
        if(windowSize < 4 || windowSize > 9){
            cout << "Wielkosc okna musi miec wartosc pomiedzy 4 a 9." << endl;
        }else{
        tmp = 1;
        vector<pair<int, pair<int, string>>> windows = getWindows(trimmedSequences, trimmedPositions, windowSize);
        }
    }
    createGraph();
    find_star(graph, windows);
}

