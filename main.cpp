
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <math.h>
#include <time.h>


template<typename T>
struct NMatrix
{
    int rowSize;
    int colSize;
    std::vector<std::vector<T>> Matrix;
    
    NMatrix(int rowSize, int colSize)
    {
        this->rowSize = rowSize;
        this->colSize = colSize;
    }
    
    int  getRowSize()
    {
        return rowSize;
    }
    
    int  getColSize()
    {
        return colSize;
    }
};

void fileReading(std::vector<double>& v)
{
    {
        std::ifstream file;
        std::string filename;
        double number;
        std::cout << "\nВведите путь до файла /Users/dina/Desktop/GAZPROM.txt : ";
        std::cin >> filename;
        file.open(filename);
        if(file.is_open())
        {
            while(file >> number){
                v.push_back(number);
            }
        }
        else
        {
            std::cout << "\nФайл " << filename << " отсутствует или нет возможности его открыть!" ;
        }
        file.close();
    }
}

void printV(std::vector<double> v)
{
    std::cout << "\n" << std::endl;

    for(int i = 0; i < v.size(); i++ )
    {
        std::cout << v.at(i) << " " ;
    }
    std::cout << "\n" << std::endl;

}

void printVN(std::vector<double> v)
{
    std::cout << "\n" << std::endl;

    for(int i = 0; i < v.size(); i++ )
    {
        std::cout << i + 1 << "\t" << v.at(i) << std::endl;
    }
    std::cout << "\n" << std::endl;

}

void preProcess(std::vector<double> &v)
{
    for(int i = 0; i < v.size(); i++)
    {
        if (i + 1 == v.size())
            break;
        v.at(i) = v.at(i + 1) - v.at(i);
    }
    v.pop_back();
}

void welcome()
{
    char ch;
    std::cout << "Добро пожаловать в нейросеть, которая предскажет вам котировку интересующих вас акций!\n" << std::endl;
    std::cout << "Для продолжения нажмите Enter" << std::endl;
    ch = std::cin.get();
    if (ch == '\n')
    {
        
    }
}

std::vector<double> getLayer(std::vector<double> full, int begin, int end)
{
    std::vector<double> v;
    for(int i = begin; i < end; i++)
    {
        v.push_back(full.at(i));
    }
    return v;
    
}

void weighFill(NMatrix<double> &W, bool pt)
{
    int rowSize = W.getRowSize();
    int colSize = W.getColSize();
    W.Matrix.resize(rowSize);
    if (pt) {
        for (int i = 0; i < rowSize; i++)
        {
            W.Matrix[i].resize(colSize);
            for(int j = 0; j < colSize; j++)
            {
                W.Matrix[i][j] = (double)(rand()) / RAND_MAX * 1 - 0.5;
    //            W.Matrix[i][j] = 1;
            }
        }
    } else {
        for (int i = 0; i < rowSize; i++)
        {
            W.Matrix[i].resize(colSize);
            for(int j = 0; j < colSize; j++)
            {
                W.Matrix[i][j] = 0;
            }
        }
    }
    
}


void printM(NMatrix<double> W)
{

    int rowSize = W.getRowSize();
    int colSize = W.getColSize();
    for(int i = 0; i < rowSize; i++)
    {
        std::cout << "\n" << std::endl;
        for(int j = 0; j < colSize; j++)
        {
            std::cout << W.Matrix[i][j] << " " ;
        }
    }
    std::cout << "\n" << std::endl;

}

std::vector<double> vecMultMat(std::vector<double> &vector, NMatrix<double> &W, int bias)
{
    
    int row = W.getRowSize();
    int col = W.getColSize();
    std::vector<double> vectorNew;
    if (vector.size() == row)
    {
        double tmp = 0.0;
        for(int j = 0; j < col; j++)
        {
            for(int i = 0; i < row; i++)
            {
                tmp += W.Matrix[i][j] * vector[i] + bias;
            }
            vectorNew.push_back(tmp);
            tmp = 0.0;
        }
    }
    else{
        std::cout << "Такие матрицы нельзя перемножить, так как количество строк вектора не равно количеству строк матрицы В." << std::endl;
    }
    return vectorNew;

}

std::vector<double> matMultVec(std::vector<double> &vector, NMatrix<double> &W)
{
    
    int row = W.getRowSize();
    int col = W.getColSize();
    std::vector<double> vectorNew;
    if (vector.size() == col)
    {
        double tmp = 0.0;
        for(int i = 0; i < row; i++)
        {
            for(int j = 0; j < col; j++)
            {
                tmp += W.Matrix[i][j] * vector[j];
            }
            vectorNew.push_back(tmp);
            tmp = 0.0;
        }
    }
    else{
        std::cout << "Такие матрицы нельзя перемножить, так как количество строк вектора не равно количеству строк матрицы В." << std::endl;
    }
    return vectorNew;

}

void activate(std::vector<double> &v)
{
    double tmp;
    for(int i = 0; i < v.size(); i++)
    {
        tmp = tanh(v[i]);
        v[i] = tmp;
        tmp = 0.0;
    }
}

void reWeight(NMatrix<double> &WM, NMatrix<double> &W, std::vector<double> delta, std::vector<double> layer1, std::vector<double> layer2, double trainSpeed, double moment)
{
    for (int i = 0; i < W.getRowSize(); i++)
    {
        for(int j = 0; j < W.getColSize(); j++)
        {
            W.Matrix[i][j] = W.Matrix[i][j] + trainSpeed * delta[j] * (1 - layer2[j] * layer2[j]) * layer1[i] + moment * WM.Matrix[i][j];
            WM.Matrix[i][j] = trainSpeed * delta[j] * (1 - layer2[j] * layer2[j]) * layer1[i] + moment * WM.Matrix[i][j];
        }
    }
}


double errorCalc(std::vector<double> idealV, int pos, std::vector<double> actualV)
{
    double error;
    double ideal = tanh(idealV.at(pos));
    return error = (ideal - actualV.at(0)) * (1 - actualV.at(0) * actualV.at(0));
}


int main(int argc, const char * argv[]) {
    
    welcome();
   

//-------Считывание данных из файла и их хранение и обработка--------------//
    std::vector<double> fullData;
    fileReading(fullData);
    std::vector<double> fullPostData(fullData);
    preProcess(fullPostData);
    
//-------Создаем слои нейросети--------------------------------------------//
    int sizeInputLayer;
    std::cout << "\nВведите количество нейронов входного слоя = ";
    std::cin >> sizeInputLayer;
    int sizeLayer1 = sizeInputLayer * 2;
    int sizeLayer2 = sizeInputLayer;
    int sizeLayer3 = sizeInputLayer / 4;
    int sizeLayer4 = sizeInputLayer / 32;
    double bias = 0.5;
    
//-------Создаем нулевой слой нейросети------------------------------------//
    std::vector<double> inputLayer;
    
//-------Создаем матрицы весов и матрицы ошибок нейросети------------------//
    NMatrix<double> W1 (sizeInputLayer, sizeLayer1);
    NMatrix<double> W2 (sizeLayer1, sizeLayer2);
    NMatrix<double> W3 (sizeLayer2, sizeLayer3);
    NMatrix<double> W4 (sizeLayer3, sizeLayer4);
    
    NMatrix<double> WM1 (sizeInputLayer, sizeLayer1);
    NMatrix<double> WM2 (sizeLayer1, sizeLayer2);
    NMatrix<double> WM3 (sizeLayer2, sizeLayer3);
    NMatrix<double> WM4 (sizeLayer3, sizeLayer4);
    
    

//-------Создаем последующие слои нейросети---------------------------------//
    std::vector<double> layer1, layer2, layer3, layer4;
    
    
//-------МЕТОД ОБРАТНОГО РАСПРОСТРАНЕНИЯ ОШИБКИ----------------------------//
    int maxEpoch;
    std::cout << "\nВведите количество эпох обучения нейронной сети = ";
    std::cin >> maxEpoch;
    std::cout << "\n";
    
    double moment;
    std::cout << "Введите постоянную момента от (-1; +1) = ";
    std::cin >> moment;
    std::cout << "\n";

    int trainSet = floor(fullPostData.size() / sizeInputLayer);
    
//-------Векторы ошибок и подсчет------------------------------------------//
    std::vector<double> layer1Error(sizeLayer1);
    std::vector<double> layer2Error(sizeLayer2);
    std::vector<double> layer3Error(sizeLayer3);
    std::vector<double> layer4Error(sizeLayer4);
    double error;
  
    
    for (double trainSpeed = 0.1; trainSpeed < 1 ; trainSpeed = trainSpeed + 0.1) {
//        double trainSpeed = 0.5;
//-------Заполняем нулевой слой нейросети------------------------------------//
        inputLayer = getLayer(fullPostData, 0, sizeInputLayer);

//-------Заполняем матрицы весов и матрицы ошибок нейросети------------------//
        weighFill(W1, true);
        weighFill(W2, true);
        weighFill(W3, true);
        weighFill(W4, true);
        
        weighFill(WM1, false);
        weighFill(WM2, false);
        weighFill(WM3, false);
        weighFill(WM4, false);
        
//-------Заполняем и  последующие слои нейросети----------------------------//
        layer1 = vecMultMat(inputLayer, W1, bias);
        activate(layer1);
        layer2 = vecMultMat(layer1, W2, bias);
        activate(layer2);
        layer3 = vecMultMat(layer2, W3, bias);
        activate(layer3);
        layer4 = vecMultMat(layer3, W4, bias);
        activate(layer4);
        error = errorCalc(fullPostData, sizeInputLayer - 1, layer4);
        for (int i = 0; i < maxEpoch; i++) {
            for (int j = 0; j < trainSet; j++) {
                layer4Error.at(0) = error;
                layer3Error = matMultVec(layer4Error, W4);
                layer2Error = matMultVec(layer3Error, W3);
                layer1Error = matMultVec(layer2Error, W2);
                reWeight(WM1, W1, layer1Error, inputLayer, layer1, trainSpeed, moment);
                reWeight(WM2, W2, layer2Error, layer1, layer2, trainSpeed, moment);
                reWeight(WM3, W3, layer3Error, layer2, layer3, trainSpeed, moment);
                reWeight(WM4, W4, layer4Error, layer3, layer4, trainSpeed, moment);

                error = errorCalc(fullPostData, (sizeInputLayer - 1) * j, layer4);

                inputLayer = getLayer(fullPostData, sizeInputLayer * j, sizeInputLayer * (j + 1));
                layer1 = vecMultMat(inputLayer, W1, bias);
                activate(layer1);
                layer2 = vecMultMat(layer1, W2, bias);
                activate(layer2);
                layer3 = vecMultMat(layer2, W3, bias);
                activate(layer3);
                layer4 = vecMultMat(layer3, W4, bias);
                activate(layer4);
            }
            std::cout << "Эпоха " << i + 1 << "/" << maxEpoch << " завершилась" << std::endl;
        }
        std::cout << "\n//---------------------------//\n" << std::endl;
        std::cout << "Для скорости обучения : " << trainSpeed << std::endl;
        std::cout << "Рассчет ошибки : " << error << std::endl;
        std::cout << "Реальная цена : 205.13" << std::endl;
        std::cout << "Прогнозируемая цена составит : " << fullData.at(fullPostData.size()) + layer4.back() << std::endl;
        std::cout << "\n";
    }
    
    std::cout << "\n";
    return 0;
}
