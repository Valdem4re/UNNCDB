/**
 * @file cif.hpp
 * @author Nikolay V. Somov (somov@phys.unn.ru)
 * @brief Библиотека для работы с CIF-файлами
 * @version 1.0
 * @date 2025-03-13
 * 
 * Данная библиотека состоит из одного заголовочного файла и не требует установки. При реализации использовались только 
 * стандартные библиотеки C++ (C++17).
 * Компонентны библиотеки:
 * CifItem - класс базового элемента. Храни в себе тег и значение.
 * CifItemLoop - класс таблицы значений. Позволяет хранить в себе множества данных. Например, таблицу атомов, таблицу длин связей и. т.д.
 * CifItemLoopFill - вспомогательный класс. Служит для удобного заполнения таблиц.
 * CifBlock - класс блока в CIF-файле. Содержит списки элементов CifItem и таблиц CifItemLoop.
 * CifBuffer - вспомогательный класс буфера данных. Служит для парсинга CIF-файла.
 * CifFile - класс CIF-файла. Позволяет читать CIF-файл и записывать. Содержит список блоков (CifBlock).
 * @copyright Copyright (c) 2025
 * 
 */
#ifndef CIF_HPP
#define CIF_HPP

//==================================================
#define CIF_PARSER_VERSION "1.0"
//==================================================

#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <exception>
#include <type_traits>
#include <cmath>

using namespace std;

//==================================================

namespace cif {

//--------------------------------------------------
enum tCifDataType
{
    CIF_NONE,
    CIF_WORD,
    CIF_STRING,
    CIF_TEXT,
    CIF_TABLE
};

//! Максимальное число "слов" в буфере чтения
#define CIF_MAX_WORD_IN_BUFFER 512

//--------------------------------------------------
inline string strip_sp(const string& s);
//--------------------------------------------------
/**
 * @brief Класс базового элемента CIF-файла. Храни в себе пару тег и значение.
 * 
 */
class CifItem
{
    public:
        tCifDataType type;
        string tag;
        string value;

        CifItem():type(CIF_NONE){

        };

        CifItem(const CifItem& item):
            type(item.type), tag(item.tag), value(item.value){

        };

        CifItem(const string& tag, const string& value)
        {
            this->tag = tag;
            this->value = value;
            auto_type();
        };

        CifItem(const string& tag, const string& value, const tCifDataType type)
        {
            this->tag = tag;
            this->value = value;
            this->type = type;
        };

        //! True если содержится пробел внутри значения (value)
        inline bool has_space() {
            for (auto c : value)
                if (c == ' ' || c == '\t') return true;
            return false;
        };

        //! True если содержит переходы на новую строку
        inline bool has_lines() {
            for (auto c : value)
                if (c == '\n') return true;
            return false;
        };

        //! Автоматическое определение типов
        tCifDataType auto_type()
        {
            type = CIF_TEXT;
            if (!has_lines())
            {
                if (has_space())
                {
                    if (value.size() < 50) type = CIF_STRING;
                }
                else
                {
                    type = (value.size()<=12)?(CIF_WORD):(CIF_STRING);
                }
            }
            return type;
        };

        CifItem& operator=(const CifItem& item)
        {
            type = item.type;
            tag = item.tag;
            value = item.value;

            return *this;
        };

        bool operator==(const CifItem& item) const
        {
          return tag == item.tag;
        };

        //! Возвращает значение в виде строки с удаленным значением в круглых скобках. "0.1234(5)" -> "0.1234"
        string as_strip_bk() const
        {
            size_t pos = value.find("(");
            if ( pos != string::npos)
            {
                return value.substr(0, pos);
            }
            else
                return value;
        };

        //! Возвращает значение как целое число
        int as_int() const
        {
            stringstream ss(as_strip_bk());
            int result;
            ss>>result;
            return result;
        };

        //! Возвращает значение как действительное число
        double as_double() const
        {
            stringstream ss(as_strip_bk());
            double result = 0.0;
            ss >> result;
            return result;
        }

        //! Сравнение тегов
        inline bool check(const string& tag) const
        {
            return this->tag == tag;
        }

        friend ostream& operator<<(ostream& out, const CifItem& item)
        {
            switch (item.type) {
            case CIF_WORD:
                out<<item.tag<<"   "<<item.value;
                break;
            case CIF_STRING:
                out<<item.tag<<"   \'"<<item.value<<"\'";
                break;
            case CIF_TEXT:
                out<<item.tag<<"\n;\n"<<item.value<<"\n;";
                break;
            default:
                out<<item.tag<<"   "<<item.value;
                break;
            }

            return out;
        };
};

//--------------------------------------------------

/**
 * @brief Класс таблицы данных (loop__)
 * 
 */
class CifItemLoop
{
    public:
        vector<string> header;
        vector<vector<string>> data;

        CifItemLoop(){

        };

        CifItemLoop(const CifItemLoop& item):
            header(item.header), data(item.data)
        {

        };

        void clear()
        {
            header.clear();
            data.clear();
        };

        string& operator()(const int i, const int j)
        {
            return data[i][j];
        };

        //! Добавляет строку значений в таблицу
        void add_row(const vector<string>& row)
        {
            data.push_back(row);
        };

        //! Добавляет строку значений в таблицу
        vector<string>& add_empty_row()
        {
            vector<string> row(header.size());
            data.push_back(row);
            return data.back();
        };

        //! Количество строк в таблице
        int rows() const
        {
            return data.size();
        };

        //! Количество столбцов в таблице
        int cols() const
        {
            if (!data.empty()) return data[0].size();
            else return 0;
        };

        //! Ссылка на строку таблицы
        vector<string>& get_row(const int index)
        {
            return data[index];
        };

        //! Ссылка на строку таблицы
        const vector<string>& get_row(const int index) const
        {
            return data[index];
        };

        CifItemLoop& operator=(const CifItemLoop& item)
        {
            header = item.header;
            data = item.data;
            return *this;
        };

        //! Проверка тега
        bool check_tag(const string& tag) const
        {
            for (auto t : header)
                if (t == tag) return true;
            return false;
        };

        //! Возвращает индекс тега в шапке таблицы (номер столбца)
        int get_header_index(const string& tag) const
        {
            const int sz = header.size();
            for (int i=0; i<sz; i++)
                if (header[i] == tag) return i;
            return -1;
        };

        friend ostream& operator<<(ostream& out, const CifItemLoop& item)
        {
            out<<"loop_"<<endl;

            for (auto h : item.header) out<<h<<endl;

            for (auto row : item.data)
            {
                out<<endl;
                for (auto v: row)
                {
                    CifItem d("", v);
                    out<<d;
                }
            }

            return out;
        };
};

//--------------------------------------------------

/**
 * @brief Вспомогательный класс. Служит для более удобного заполнения таблицы
 * 
 */
class CifItemLoopFill
{
    protected:

        vector<string> row;
        CifItemLoop& loop;

    public:

        CifItemLoopFill(CifItemLoop& _loop):loop(_loop)
        {
        };

        ~CifItemLoopFill()
        {
        };

        void add(const string& value)
        {
            row.push_back(value);

            if (row.size() == loop.header.size())
            {
                loop.data.push_back(row);
                row.clear();
            }
        };

        CifItemLoopFill& operator<<(const string& value)
        {
            add(value);
            return *this;
        }

        bool check_loop()
        {
            return row.empty();
        };

        void check_loop_throw()
        {
            if (!row.empty())
            {
                stringstream ss;
                ss<<loop;
                throw runtime_error("CIF ERROR: Not complete data set for loop. " + ss.str());
            }
        };
};

//--------------------------------------------------

/**
 * @brief Класс блока данных в CIF-файле.
 * Содержит списки элементов и таблиц значений.
 */

class CifBlock
{
    public:
        vector<CifItem> items;
        vector<CifItemLoop> loops;
        string name;

        CifBlock()
        {

        };

        CifBlock(const string& name):name(name)
        {

        };

        CifBlock(const CifBlock& block):
            items(block.items), loops(block.loops), name(block.name)
        {

        };

        void clear()
        {
            name="";
            items.clear();
            loops.clear();
        };

        //! Добавить элемент
        void add_item(const CifItem& item)
        {
            items.push_back(item);
        };

        //! Добавить элемент
        void add_item(const string& tag, const string& value)
        {
            CifItem item(tag, value);
            items.push_back(item);
        };

        //! Добавить таблицу
        void add_loop(const CifItemLoop& loop)
        {
            loops.push_back(loop);
        };


        //! найти элемент по тегу
        CifItem* find_item(const string& tag)
        {
            for (CifItem& item : items)
                if (item.check(tag)) return &item;
            return 0;
        };

        //! найти элемент по тегу
        const CifItem* find_item(const string& tag) const
        {
            for (const CifItem& item : items)
                if (item.check(tag)) return &item;
            return 0;
        };

        //! найти таблицу по тегу
        CifItemLoop* find_loop(const string& tag)
        {
            for (CifItemLoop& item : loops)
                if (item.check_tag(tag)) return &item;
            return 0;
        };

        //! найти таблицу по тегу
        const CifItemLoop* find_loop(const string& tag) const
        {
            for (const CifItemLoop& item : loops)
                if (item.check_tag(tag)) return &item;
            return 0;
        };

        friend ostream& operator<<(ostream& out, const CifBlock& block)
        {
            out<<"data_"<<block.name<<endl<<endl;

            for (auto x : block.items) out<<x<<endl;
            for (auto x : block.loops) out<<x<<endl;

            return out;
        };
};

//--------------------------------------------------

/**
 * @brief Вспомогательный класс буфера данных.
 * Выполняет последовательное чтение данных из потока (CIF-файла) с
 * последующим разбором (парсингом).
 */

class CifBuffer: public queue<string>
{
    private:
        bool __remove;
        bool __eof;
        istream* lpFile;
    public:

        CifBuffer():__remove(false), __eof(true), lpFile(0)
        {

        };

        CifBuffer(const string& path)
        {
            ifstream* __f = new ifstream(path);
            lpFile = __f;
            
            __remove = true;

            if (lpFile->bad())
            {
                __eof = true;
            }
            else
            {
                __eof = false;
            }

            if (!__f->is_open())
            {
                stringstream ss;
                ss<<"Не могу открыть файл: "<<path<<". Возможно файл не существует.";
                throw runtime_error(ss.str());
            }
        };

        CifBuffer(istream& input)
        {
            lpFile = &input;
            __remove = false;
            __eof = false;
        };

        ~CifBuffer()
        {
            if (lpFile && __remove) delete lpFile;
        };

        bool open(istream& input)
        {
            close();

            lpFile = &input;
            __remove = false;
            __eof = false;
            return true;
        };

        bool open(const string& path)
        {
            close();

            lpFile = new ifstream(path);
            __remove = true;
            __eof = false;

            return (bool)(*lpFile);
        };

        void close()
        {
            if (lpFile && __remove) delete lpFile;
            lpFile = 0;
            __remove = false;
        };

        bool get(string& s)
        {
            if (empty())
            {
                if (_read_data())
                    return get(s);
                else
                    return false;
            }
            else
            {
                return _read_value(s);
            }
        };

        void pop()
        {
            queue<string>::pop();
        };

        bool empty()
        {
            return queue<string>::empty();
        };

    protected:

        bool _read_data()
        {
            int count = CIF_MAX_WORD_IN_BUFFER;
            string buf;

            while ((!__eof) && (count>0))
            {
                std::getline(*lpFile, buf);
                
                if (buf.empty() || buf[0] == '#') 
                {
                    __eof = lpFile->eof();
                    continue;
                }
                
                if (buf[0] =='_' && buf[buf.size()-1] == ';')
                {
                    buf[buf.size()-1] = '-';
                }

                stringstream ss(buf);
                string s;

                if (buf.size()>0 && buf[0] == ';')
                {
                    string text = strip_sp(buf.substr(1));
                    do {
                        std::getline(*lpFile, buf);
                        __eof = lpFile->eof();
                        if (buf[0] == ';') break;
                        else 
                        {
                            if (text.empty()) text += buf;
                            else
                                text += "\n" + buf;
                        }
                    } while (!__eof);

                    this->push(text);
                    count--;

                    continue;
                }
                
                while (ss >> s)
                {
                    if (!s.empty())
                    {
                        this->push(s);
                        count--;
                    }
                }

                __eof = lpFile->eof();
            }

            return (!__eof) || (!empty());
        };

        inline bool _check_semicolon(const string&s, string& addon_)
        {
            if (s.size() > 0 && s[0] == ';')
            {
                addon_ = cif::strip_sp(s.substr(1));
                return true;
            }
            else 
                return false;
        };

        inline bool _read_value(string& s)
        {
            s = front();
           
            if (s == ";")
            {
                pop();
                string t;
                s.clear();

                while ( !(empty() && __eof) )
                {
                    t = front();

                    if (t == ";") break;
                    else
                    {
                        if (!s.empty()) s+=" ";
                        s+=t;
                    }

                    pop();

                    if (empty()) _read_data();
                }

                front() = s;
            }
            else if (s[0] == '\'' || s[0] == '\"')
            {
                char bk = s[0];
                s.erase(s.begin());

                if (s.size() > 0 && s[s.size()-1] == bk)
                {
                    s.pop_back();
                    return true;
                }
                
                pop();

                string t;

                while ( !(empty() && __eof) )
                {
                    if (empty())  
                    {
                        _read_data(); 
                        continue;
                    }
                
                    t = front();

                    if (t[t.size()-1] == bk)
                    {
                        t.pop_back();
                        if (!s.empty()) s+=" ";
                        s+=t;
                        break;
                    }
                    else
                    {
                        if (!s.empty()) s+=" ";
                        s+=t;
                    }

                    pop();
                    if (empty()) _read_data();
                }
                front() = s;
            }

            return true;
        };
};

//--------------------------------------------------
/**
 * @brief Класс CIF-файла. Содержит список блоков.
 * Поддерживает режим полного чтения CIF-файла и режим последовательного (поблочного) чтения
 * для обработки больших файлов
 */
class CifFile: public vector<CifBlock>
{
    protected:

        CifBuffer _buf;

    public:

        CifFile()
        {

        };

        CifFile(const CifFile& cif):
            vector<CifBlock>(cif), _buf(cif._buf)
        {

        };

        CifFile(const string& path):
            vector<CifBlock>(), _buf(path)
        {

        };

        CifFile(istream& input):
            vector<CifBlock>(), _buf(input)
        {

        };

        ~CifFile(){};

        /**
          Низкоуровневая функция чтения одного блока из файла.
          Применяется для чтения длинных файлов (последовательный доступ).
          При успешном выполнении чтения возвращает true.
        */
        bool read_block(CifBlock& block)
        {
            block.clear();
            string w;

            while (_buf.get(w))
            {
                if (w.find("data_") == 0)
                {
                    if (block.name.empty())
                    {
                        block.name = w.substr(5, w.size() - 5);
                        _buf.pop();
                    }
                    else break;
                    continue;
                }

                if (w[0] == '_')
                {
                    _buf.pop();
                    string tag = w;

                    if (_buf.get(w))
                    {
                        string value = w;
                        block.add_item(tag, value);
                        _buf.pop();
                    }
                    else
                    {
                        throw runtime_error("CIF ERROR: No value found for tag ["+tag+"].");
                    }
                    continue;
                }

                if (w == "loop_") 
                {
                    _read_loop(block);
                    continue;
                }

                throw runtime_error("Invalid CIF-file format: " + w + ".");
            }

            return block.items.size() > 0;
        };

        /**
            Читает все данные из cif-файла. Возвращает число прочитанных блоков.
        */
        int read_blocks()
        {
            CifBlock block;

            while (read_block(block))
            {
                this->push_back(block);
            }

            return this->size();
        };

    private:

        inline bool _read_loop(CifBlock& block)
        {
            CifItemLoop loop;

            _buf.pop();

            // read header
            string w;

            while (_buf.get(w))
            {
                if (w[0] == '_')
                {
                    loop.header.push_back(w);
                    _buf.pop();
                }
                else break;
            }

            if (loop.header.empty())
                throw runtime_error("CIF ERROR: Loop is empty. Line: " + w + "...");

            // read data

            CifItemLoopFill wr(loop);

            while (_buf.get(w))
            {
                if (w[0] == '_' || w == "loop_" || w.substr(0, 5)=="data_") break;
                else
                {
                    wr<<w;
                    _buf.pop();
                }
            }

            wr.check_loop_throw();

            block.add_loop(loop);

            return wr.check_loop();
        }
};

//--------------------------------------------------
//! Remove esd value in ()
inline string strip_bk(const string& value)
{
    size_t pos = value.find("(");
    if ( pos != string::npos)
    {
        return value.substr(0, pos);
    }
    else
        return value;
};

//--------------------------------------------------
//! Convert value with esd to numerical
template<typename T> T str_to_value(const string& value)
{
    stringstream ss(strip_bk(value));
    T result;
    ss>>result;
    return result;
};

//--------------------------------------------------

/**
 * @brief Возвращает значение по тегу.
 * 
 * @param block блок данных
 * @param tag тег для поиска
 * @param default_value значение, возвращаемое если тег не найден
 * @return результат поиска
 */
template<typename T> T get_value(const CifBlock& block, const string& tag, const T& default_value)
{
    const CifItem* item = block.find_item(tag);
    if (item)
    {
        return str_to_value<T>(item->value);
    }
    else
        return default_value;
};
//--------------------------------------------------

/**
 * @brief Возвращает значение по тегу.
 * 
 * @param block блок данных
 * @param tag тег для поиска
 * @param default_value значение, возвращаемое если тег не найден
 * @return string результат поиска
 */
inline string get_value(const CifBlock& block, const string& tag, const string& default_value)
{
    const CifItem* item = block.find_item(tag);
    if (item)
    {
        return item->value;
    }
    else
        return default_value;
};
//--------------------------------------------------

//! Удаляет пробельные символы в начале и конце строки
inline string strip_sp(const string& s)
{
    int first = 0;
    int sz = s.size();
    int last = sz - 1;
    while (first < sz && s[first] <= ' ') first++;
    while (last>=0 && s[last]  <= ' ') last--;

    string r = s.substr(first, last);

    return r;
}
//--------------------------------------------------
/**
 * @brief Extract value and esd
 * 
 * @param s 
 * @param value 
 * @param esd 
 * @return true 
 * @return false 
 */
inline bool extract_values(const string& s, double& value, double& esd)
{
    value = esd = 0;
    int state = 0;
    int p_point = -1;
    int p_first = -1;
    int p_bk1 = -1;
    int p_bk2 = -1;
    size_t sz = s.size();

    for (size_t i=0; i < sz; ++i)
    {
        if ( p_first < 0 && s[i] > ' ' ) p_first = i;

        switch (s[i])
        {
            case '.':
                state = 1;
                p_point = i;
                break;
            case '(':
                state = 2;
                p_bk1 = i;
                break;
            case ')':
                state = 3;
                p_bk2 = i;
                break;
        }
    }

    if (p_first >= 0)
    {
        size_t n_ = (state>=2) ? (p_bk1 - p_first) : (sz);
        stringstream sv(s.substr(p_first, n_));
        sv >> value;
    }
    else 
        return false;

    if (p_bk1>0 && p_bk2>0) 
    {
        size_t n_ = p_bk2 - p_bk1 - 1;
        stringstream se(s.substr(p_bk1+1, n_));
        se >> esd;
        esd = esd * std::pow(10, -(p_bk1 - p_point - 1));
    }
    
    return true;
}

}; // namespace cif

#endif // CIF_HPP
