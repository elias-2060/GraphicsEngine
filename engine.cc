#define _USE_MATH_DEFINES
#include "easy_image.h"
#include "ini_configuration.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include "cmath"
#include "l_parser.h"
#include "stack"

using namespace img;
using namespace std;

EasyImage draw2DLines(Lines2D &lines, const int size, Color &backgroundColor){
    double xMin = lines.begin()->p1.x;
    double xMax = lines.begin()->p1.x;
    double yMin = lines.begin()->p1.y;
    double yMax = lines.begin()->p1.y;

    /// kleinste en grootste waarde voor x en y zoeken tussen alle punten
    for (const auto& line : lines){
        if (line.p1.x < xMin)
            xMin = line.p1.x;
        if (line.p1.x > xMax)
            xMax = line.p1.x;
        if (line.p1.y < yMin)
            yMin = line.p1.y;
        if (line.p1.y > yMax)
            yMax = line.p1.y;

        if (line.p2.x < xMin)
            xMin = line.p2.x;
        if (line.p2.x > xMax)
            xMax = line.p2.x;
        if (line.p2.y < yMin)
            yMin = line.p2.y;
        if (line.p2.y > yMax)
            yMax = line.p2.y;
    }
    double xRange = xMax - xMin;
    double yRange = yMax - yMin;
    double imageX = size * (xRange / (max(xRange, yRange)));
    double imageY = size * (yRange / (max(xRange, yRange)));
    double d = 0.95 * (imageX / xRange);
    EasyImage image(imageX,imageY, backgroundColor);

    for (auto &line : lines) {
        line.p1.x *= d;
        line.p1.y *= d;
        line.p2.x *= d;
        line.p2.y *= d;
    }

    double DCx = d * ((xMin + xMax) / 2);
    double DCy = d * ((yMin + yMax) / 2);
    double dX = imageX / 2 - DCx;
    double dY = imageY / 2 - DCy;

    for (auto &line : lines) {
        line.p1.x += dX;
        line.p1.y += dY;
        line.p2.x += dX;
        line.p2.y += dY;
    }
    /// draw the lines
    for (const auto& line : lines){
        image.draw_line(::lround(line.p1.x), ::lround(line.p1.y), ::lround(line.p2.x), ::lround(line.p2.y), Color(
                line.color.red, line.color.green, line.color.blue));

    }
    return image;
}

string getReplacements(const string &str, const LParser::LSystem2D &parser, int iterations){
    if (iterations == 0)
        return str;

    string newString;

    for (char j: str) {
        if (j == '+' or j == '-' or j == '(' or j == ')')
            newString += j;
        else
            newString += parser.get_replacement(j);
    }
    iterations--;
    return getReplacements(newString, parser, iterations);
}

EasyImage draw2DLSystem(const LParser::LSystem2D &parser, const int size, const vector<double>& lineColor, const vector<double>& backgroundcolor){
    Lines2D lines;
    string strings = getReplacements(parser.get_initiator(), parser, parser.get_nr_iterations());
    double angle = parser.get_starting_angle() * (M_PI / 180);
    Point2D pos(0,0);
    stack<vector<double>> stack;
    for (char j: strings){
        if (j == '+'){
            angle += parser.get_angle() * (M_PI / 180);
        }
        else if (j == '-'){
            angle -= parser.get_angle() * (M_PI / 180);
        }
        else if (j == '('){
            vector<double> temp;
            temp.push_back(pos.x);
            temp.push_back(pos.y);
            temp.push_back(angle);
            stack.push(temp);
        }
        else if (j == ')'){
            vector<double> temp = stack.top();
            stack.pop();
            pos.x = temp[0];
            pos.y = temp[1];
            angle = temp[2];
        }
        else{
            Point2D prevPos = pos;
            pos.x += cos(angle);
            pos.y += sin(angle);
            if (parser.draw(j)){
                Color color(lineColor[0], lineColor[1], lineColor[2]);
                lines.emplace_back(prevPos, pos, Color(color));
            }
        }
    }

    Color c(backgroundcolor[0], backgroundcolor[1], backgroundcolor[2]);
    return draw2DLines(lines, size, c);
}


img::EasyImage generate_image(const ini::Configuration &configuration){
    EasyImage image;
    string type = configuration["General"]["type"].as_string_or_die();
    int size = configuration["General"]["size"].as_int_or_die();
    vector<double> backgroundcolor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();

    if (type == "2DLSystem") {
        string inputfile = "files/";
        inputfile += configuration["2DLSystem"]["inputfile"].as_string_or_die();
        vector<double> lineColor = configuration["2DLSystem"]["color"].as_double_tuple_or_die();
        LParser::LSystem2D parser;
        ifstream input_stream(inputfile);
        input_stream >> parser;
        input_stream.close();
        image = draw2DLSystem(parser, size, lineColor, backgroundcolor);
    }
	return image;
}

int main(int argc, char const* argv[])
{
        int retVal = 0;
        try
        {
                std::vector<std::string> args = std::vector<std::string>(argv+1, argv+argc);
                if (args.empty()) {
                        std::ifstream fileIn("filelist");
                        std::string filelistName;
                        while (std::getline(fileIn, filelistName)) {
                                args.push_back(filelistName);
                        }
                }
                for(std::string fileName : args)
                {
                        ini::Configuration conf;
                        try
                        {
                                std::ifstream fin(fileName);
                                if (fin.peek() == std::istream::traits_type::eof()) {
                                    std::cout << "Ini file appears empty. Does '" <<
                                    fileName << "' exist?" << std::endl;
                                    continue;
                                }
                                fin >> conf;
                                fin.close();
                        }
                        catch(ini::ParseException& ex)
                        {
                                std::cerr << "Error parsing file: " << fileName << ": " << ex.what() << std::endl;
                                retVal = 1;
                                continue;
                        }

                        img::EasyImage image = generate_image(conf);
                        if(image.get_height() > 0 && image.get_width() > 0)
                        {
                                std::string::size_type pos = fileName.rfind('.');
                                if(pos == std::string::npos)
                                {
                                        //filename does not contain a '.' --> append a '.bmp' suffix
                                        fileName += ".bmp";
                                }
                                else
                                {
                                        fileName = fileName.substr(0,pos) + ".bmp";
                                }
                                try
                                {
                                        std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                                        f_out << image;

                                }
                                catch(std::exception& ex)
                                {
                                        std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                                        retVal = 1;
                                }
                        }
                        else
                        {
                                std::cout << "Could not generate image for " << fileName << std::endl;
                        }
                }
        }
        catch(const std::bad_alloc &exception)
        {
    		//When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
    		//Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
    		//(Unless of course you are already consuming the maximum allowed amount of memory)
    		//If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
		//mark the test as failed while in reality it just needed a bit more memory
                std::cerr << "Error: insufficient memory" << std::endl;
                retVal = 100;
        }
        return retVal;
}
