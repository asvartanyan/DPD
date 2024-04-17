#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <thread>
#include<chrono>
#include<random>
#include<time.h>
using namespace std;

struct Bonds
{
 public:
 int left;
 int right;
};
class PairsOfAtoms
{
 public:
 int left;
 int right;
 bool state_of_b;
};

class Atom
{
 public: 
 vector<float> F;
 vector<float> V;
 vector<float> Vstr;
 vector<float> Fstr;

 vector<float> position;//координаты атома
 Atom()  
  {
   position.resize(3,0);
   F.resize(3,0);
   V.resize(3,0);
   Vstr.resize(3,0);
   Fstr.resize(3,0);
  }
};

void ClearOutputFile(string path_file)
{
  ofstream clr_file;
  clr_file.open(path_file,ios_base::trunc);
  clr_file.close();
}

class Denrimers
{
 public:
  string num_dend = ""; //переменная номера   дендримера
  vector<Atom> atms;  //вектор атомов (количество атомов зависит от конкретного дендримера)
  vector<Bonds> p_bonds;
  vector<PairsOfAtoms> p_atoms;
  public: 
  Denrimers(string n)       //конструктор класса
{
    num_dend = n;
};
//метод чтения координат из файла COORD
 void read_cb() 
 {
   string file_path_c = "Denrimers/G=";
   string file_path_b = file_path_c;
   file_path_c += num_dend;
   file_path_c +="/COORD";
   file_path_b +=num_dend;
   file_path_b +="/BONDS";
   
   ifstream in(file_path_c);
   ifstream inb(file_path_b);
   cout<<file_path_c;
   if(in.is_open())
   {
     cout<<"\nОткрыт COORD "<<num_dend<<endl;
     string tmp;
     int num_atoms;
     in>>tmp>>num_atoms;
     //cout<<num_atoms;
     atms.resize(num_atoms);
     //считывание координат
     for(int i = 0; i<num_atoms;i++)
     {
       //atms[i].position.resize(3);
       in>>atms[i].position[0]>>atms[i].position[1]>>atms[i].position[2];
     }
     ///////// работа с файлом BONDS
     if(inb.is_open()) //открытие файла bonds
     {
      int bonds_size;
      cout<<"\nОткрыт BONDS "<<num_dend<<endl;
      inb>>tmp>>bonds_size>>tmp>>tmp;//пропуск и считывание количества bonds с первой строки
      p_bonds.resize(bonds_size);
      while(!inb.eof()) //пока не конец файла bonds
      {
       for(int i = 0; i<bonds_size;i++) //сколько связей столько и итераций
       {
         inb>>p_bonds[i].left>>p_bonds[i].right;
       }
       break;
      }
     }
     else
      {
        cout<<"\n Не открыт BONDS "<<num_dend<<endl;
      }
   }
   else
   {
     cout<<"\nНе открыт COORD " << num_dend<<endl;
   }
   in.close();
   inb.close();
 }
//моделирование
 void alg(float dt,int nstep,double aij, float Kbond,float Lbond,float Rcut,float gamma,float sigma,int print_to_pdb)
{
  random_device rd{};
  mt19937 gen{rd()};
  normal_distribution<> d{0,1};
  ofstream out_c;
  ClearOutputFile("output_coords.pdb");
//пары с повторениями
//заполнить пары атомов (исключая пару атома с самим собой(1-1)) и выставить им состояние связи

  for(int i = 0; i<atms.size();i++)
  {
    //j = i+1;
    for(int j = 1 ;j<=atms.size();j++)
    {
     PairsOfAtoms now;
     now.left = i+1;
     now.right = j;
     if(now.left==now.right)
     {
       continue;
     }
     for(int z = 0; z<p_bonds.size();z++)
     {       
       if((now.left==p_bonds[z].left && now.right == p_bonds[z].right )|| (now.left == p_bonds[z].right && now.right == p_bonds[z].left))
       {
         now.state_of_b = true;
         break;
       }
       else
       {
           now.state_of_b = false;
       }
     }
     p_atoms.push_back(now);
    }
    
  }
/*
//Пары без повторений 
  int j = 1;
        for (int i = 0; i <= atms.size(); i++)
        {
            PairsOfAtoms now;
            now.left = i + 1; //номер текущего атома
            j++;
            for (int k = j; k <=atms.size(); k++)
            {
                now.right = k;
                if (now.left == now.right)
                {
                    continue;
                }
                for (int z = 0; z < p_bonds.size(); z++)
                {
                    if (now.left == p_bonds[z].left && now.right == p_bonds[z].right)
                    {
                        now.state_of_b = true;
                        break;
                    }
                    else
                    {
                        now.state_of_b = false;
                    }
                }
              p_atoms.push_back(now);
            }
           
        }
*/
//главный цикл по числу шагов моделирования
 for(int g = 0; g<nstep;++g) 
 {
//перемещние координат атомов в дендримере и оценка скоростей
 for(int i = 0; i<atms.size();++i)
 {
  atms[i].position[0] = atms[i].position[0]+((atms[i].V[0]* dt) +((1/2.0)*(dt * dt)*atms[i].F[0]));
  atms[i].position[1] = atms[i].position[1]+((atms[i].V[1]* dt)+((1/2.0)*(dt * dt)*atms[i].F[1]));
  atms[i].position[2] = atms[i].position[2]+((atms[i].V[2]* dt)+((1/2.0)*(dt * dt)*atms[i].F[2]));
 
  //atms[i].Vstr[0] = atms[i].V[0]+((1/2.0)* dt *atms[i].F[0]);
  //atms[i].Vstr[1] = atms[i].V[1]+((1/2.0)* dt *atms[i].F[1]);
  //atms[i].Vstr[2] = atms[i].V[2]+((1/2.0)* dt *atms[i].F[2]);
 }
//Получение сил

int pa = 0;
int dubl_k = 0;
for(int k = 0; k<atms.size();++k) 
{
 float Rij = 0;
 float Xij = 0;
 float Yij = 0;
 float Zij = 0;
 float FbondRij = 0;
 float WD = 0;
 float WR = 0;
 float omega = d(gen); 
 vector<float> FB(3,0);
 vector<float> FC(3,0);
 vector<float> FD(3,0);
 vector<float> FR(3,0);
 vector<float> Vij(3,0);
vector<float> sum_FB(3,0);
 vector<float> sum_FC(3,0);
 vector<float> sum_FD(3,0);
 vector<float> sum_FR(3,0);
 float rv_scalar;
   for(int b = 0; b < atms.size()-1; ++b) //рассматриваем одну частицу с остальными
   {
    //расстояние между i и j частицей
      Xij = (atms[(p_atoms[pa].right-1)].position[0] - atms[p_atoms[pa].left-1].position[0]);
      Yij = (atms[p_atoms[pa].right-1].position[1] - atms[p_atoms[pa].left-1].position[1]);
      Zij = (atms[p_atoms[pa].right-1].position[2] - atms[p_atoms[pa].left-1].position[2]);

      Vij[0] = atms[p_atoms[pa].left-1].V[0] - atms[p_atoms[pa].right-1].V[0];
      Vij[1] = atms[p_atoms[pa].left-1].V[1] - atms[p_atoms[pa].right-1].V[1];
      Vij[2] = atms[p_atoms[pa].left-1].V[2] - atms[p_atoms[pa].right-1].V[2];

      Rij = sqrt((Xij*Xij)+(Yij*Yij)+(Zij*Zij));
    //WD
     if(Rij<Rcut)
     {
     WD = (Rcut-Rij)*(Rcut-Rij);
     } 
     else if(Rij>=Rcut)
     {
      WD=0;
     }
    //WR
     if(Rij<Rcut)
     {
      WR = abs(Rcut-Rij);
     } else if(Rij>=Rcut)
     {
      WR = 0;
     }
    //разветление связанная ли текущая пара или нет
    if(p_atoms[pa].state_of_b == true)
    {
      //FB
       FbondRij = (Kbond*(Rij-Lbond))/(Rij);
       FB[0] = FbondRij*Xij; //X
       FB[1] = FbondRij*Yij; //Y
       FB[2] = FbondRij*Zij; //Z
    }
    else if(p_atoms[pa].state_of_b == false)
    {
    //FC
     if(Rij<Rcut)
     {
      FC[0] = -(aij*(1-(Rij/Rcut)))*(Xij/Rij); //X
      FC[1] = -(aij*(1-(Rij/Rcut)))*(Yij/Rij); //Y
      FC[2] = -(aij*(1-(Rij/Rcut)))*(Zij/Rij); //Z
     }
     if(Rij>=Rcut)
     {
      FC[0] = 0; //X
      FC[1] = 0; //Y
      FC[2] = 0; //Z
     }
    //FD
    rv_scalar = Vij[0]*Xij/Rij + Vij[1]*Yij/Rij + Vij[2]*Zij/Rij;
    FD[0] = -gamma*WD*rv_scalar*Xij/Rij;
    FD[1] = -gamma*WD*rv_scalar*Yij/Rij;
    FD[2] = -gamma*WD*rv_scalar*Zij/Rij;
    //FR
     omega = d(gen);
      FR[0] = sigma*sqrt(1.0/dt)*WR*omega*Xij; //X
      FR[1] = sigma*sqrt(1.0/dt)*WR*omega*Yij; //Y
      FR[2] = sigma*sqrt(1.0/dt)*WR*omega*Zij; //Z
    }
  //SUM
    sum_FB[0] += FB[0];
    sum_FB[1] += FB[1];
    sum_FB[2] += FB[2];

    sum_FC[0] += FC[0];
    sum_FC[1] += FC[1];
    sum_FC[2] += FC[2];

    sum_FD[0] += FD[0];
    sum_FD[1] += FD[1];
    sum_FD[2] += FD[2];
   
    sum_FR[0] += FR[0];
    sum_FR[1] += FR[1];
    sum_FR[2] += FR[2];
   
    pa++;
  }
//расчет сил + получение скорости + запоминание сил для нового шага

  atms[k].Fstr[0] = sum_FB[0] + sum_FC[0] + sum_FD[0] + sum_FR[0];
  atms[k].Fstr[1] = sum_FB[1] + sum_FC[1] + sum_FD[1] + sum_FR[1];
  atms[k].Fstr[2] = sum_FB[2] + sum_FC[2] + sum_FD[2] + sum_FR[2];
if(dubl_k < 1)
{
  k = 0;
  dubl_k++;
  continue;
}
  atms[k].V[0] = (1.0/2)*(atms[k].F[0]+ atms[k].Fstr[0]) * dt;
  atms[k].V[1] =  (1.0/2)*(atms[k].F[1]+ atms[k].Fstr[1]) * dt;
  atms[k].V[2] =  (1.0/2)*(atms[k].F[2]+ atms[k].Fstr[2]) * dt;

  atms[k].F[0] = atms[k].Fstr[0];
  atms[k].F[1] = atms[k].Fstr[1];
  atms[k].F[2] = atms[k].Fstr[2];
}
//вывод в pdb

 if(g % print_to_pdb == 0) 
 {
 FILE *hfile = fopen("output_coords.pdb", "a+");
    fprintf(hfile, "COMPND    UNNAMED\n");
    fprintf(hfile, "AUTHOR    GENERATED BY TvGU\n");
    fprintf(hfile,"MODEL        0\n");
    for(int i = 0; i<atms.size();++i)
    {
    
    fprintf(hfile, "ATOM%7i    %c%4s A%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           %c\n", 
	     (i+1),'N', "ABSD", 1, atms[i].position[0], atms[i].position[1], atms[i].position[2], 1.00, 0.00, 'N');
    }
    for(int j = 0 ; j<p_bonds.size();j++)
    {
      fprintf(hfile, "CONECT%5d%5d\n", p_bonds[j].left, p_bonds[j].right);
    }
    fprintf(hfile,"TER\nENDMDL\n");
    fclose(hfile);
  }
}
}
};
void PrintPairs(Denrimers dend)
{
  for(int i = 0; i<dend.p_atoms.size();i++)
  {
    cout<<dend.p_atoms[i].left<<"\t"<<dend.p_atoms[i].right<<"\t"<<dend.p_atoms[i].state_of_b<<endl;
  }
}
int main() 
{
srand(time(NULL));

  int nstep,print_to_pdb = 0;
  double a = 1.0;
  float gamma,dt,KbT,Lbond,Kbond,Rcut = 0;
  float sigma = 0;

 Denrimers g_1 = Denrimers("1");
 Denrimers g_2 = Denrimers("2");
 Denrimers g_3 = Denrimers("3");
 Denrimers g_4 = Denrimers("4");
//чтение из файла input_file
  std::ifstream in("input_file"); 
  if(in.is_open())  //условие найденого файла
  {
    std::cout<<"\nFile opened!\n";
    string tmp;
    in>>tmp>>a>>tmp>>gamma>>tmp>>dt>>tmp>>nstep>>tmp>>KbT>>tmp>>Lbond>>tmp>>Kbond>>tmp>>Rcut>>tmp>>print_to_pdb;
   
  }
    else //услование не найденного файла
     {
     std::cout<<"\nFile is not opened!\n";
     }
  in.close();
  cout << a << "\n" << gamma << "\n" << dt << "\n" << nstep << "\n" << KbT << "\n" << Lbond << "\n" << Kbond << "\n" << Rcut << "\n" << print_to_pdb << "\n";

sigma = sqrt(2*gamma*KbT);

 clock_t start = clock();
 g_3.read_cb();
 g_3.alg(dt,nstep,a,Kbond, Lbond, Rcut, gamma, sigma, print_to_pdb);

 PrintPairs(g_3);
  clock_t end = clock();
  double seconds = (double)(end - start) / CLOCKS_PER_SEC;
  cout<<"\nВремя работы: \t"<<seconds*1000<< " m/s\n";
  return 0;
}