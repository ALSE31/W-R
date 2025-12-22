#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <limits>


const int max_age_of_rabbits = 3;
const int max_age_of_wolfs = 9;
const int size_of_iland = 120;
const int start_count_of_rabbits = 18;
const int start_count_of_wolfs = 18;
const int max_rabbits_count = 700;
const int max_wolf_count = 250;

const bool is_searching = true;
const double quality_of_search = 0.6;

const int STEPS = 1000;

struct Position {
    int x = 0;
    int y = 0;

    Position(int ix = 0, int iy = 0) : x(ix), y(iy) {};

    bool operator==(const Position& other) const noexcept {
        return (x == other.x) && (y == other.y);
    }

    bool operator!=(const Position& other) const noexcept {
        return !(*this == other);
    }

    Position& operator+=(const Position& rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    Position& operator-=(const Position& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    double distanceTo(const Position& other) const {
        int dx = x - other.x;
        int dy = y - other.y;
        return std::sqrt(dx * dx + dy * dy);
    }
};

Position operator+(Position lhs, Position rhs) {
    return lhs += rhs;
}

Position operator-(Position lhs, Position rhs) {
    return lhs -= rhs;
}

using Cord = Position;
using Point = Position;

class Animal {
public:
    Animal(int startX, int startY) : position_(startX, startY), age(0) {}
    Animal(Point startP) : position_(startP), age(0) {}
    virtual ~Animal() = default;

    virtual void Move(int width, int height) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(-1, 1);

        Position move(dis(gen), dis(gen));
        position_ += move;

        // Обеспечиваем тороидальное движение (выход за границу = появление с другой стороны)
        position_.x = (position_.x + width) % width;
        position_.y = (position_.y + height) % height;

        age++;
    }

    int getX() const { return position_.x; }
    int getY() const { return position_.y; }
    Position getPosition() const { return position_; }
    int getAge() const { return age; }
    void resetAge() { age = 0; }

protected:
    Point position_;
    int age;
};

class Rabbit : public Animal {
public:
    Rabbit(int x, int y) : Animal(x, y) {}
};

class Wolf : public Animal {
public:
    Wolf(int x, int y) : Animal(x, y) {}

    void MoveWithoutSearch(int width, int height) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(-1, 1);

        Position move(dis(gen), dis(gen));
        position_ += move;

        // Обеспечиваем тороидальное движение (выход за границу = появление с другой стороны)
        position_.x = (position_.x + width) % width;
        position_.y = (position_.y + height) % height;

        age++;
    }

    void MoveWithSearch(int width, int height, const std::vector<Rabbit>& rabbits) {
        // Если есть зайцы, пытаемся двигаться к ближайшему
        if (!rabbits.empty()) {
            // Находим ближайшего зайца
            Position target = findClosestRabbit(rabbits);

            // Вычисляем направление к цели
            int dx = 0, dy = 0;
            if (target.x > position_.x) dx = 1;
            else if (target.x < position_.x) dx = -1;

            if (target.y > position_.y) dy = 1;
            else if (target.y < position_.y) dy = -1;

            // С вероятностью двигаемся к цели, иначе случайно
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> prob(0.0, 1.0);

            if (prob(gen) < quality_of_search) {
                // Двигаемся к цели
                position_.x += dx;
                position_.y += dy;
            }
            else {
                // Случайное движение
                std::uniform_int_distribution<> dis(-1, 1);
                position_.x += dis(gen);
                position_.y += dis(gen);
            }
        }
        else {
            // Нет зайцев - случайное движение
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> dis(-1, 1);

            position_.x += dis(gen);
            position_.y += dis(gen);
        }

        // Обеспечиваем тороидальное движение
        position_.x = (position_.x + width) % width;
        position_.y = (position_.y + height) % height;

        age++;
    }

private:
    Position findClosestRabbit(const std::vector<Rabbit>& rabbits) const {
        Position closest = rabbits[0].getPosition();
        double minDistance = position_.distanceTo(closest);

        for (const auto& rabbit : rabbits) {
            double distance = position_.distanceTo(rabbit.getPosition());
            if (distance < minDistance) {
                minDistance = distance;
                closest = rabbit.getPosition();
            }
        }

        return closest;
    }
};

class Iland {
public:
    Iland(int width, int height, bool search) : n(width), m(height), search_(search) {}

    void addRabbit(int x, int y) {
        rabbits_.emplace_back(x, y);
    }

    void addWolf(int x, int y) {
        wolves_.emplace_back(x, y);
    }

    void SetSearchMod(bool search) { search_ = search; }

    void simulateStep() {
        // Движение зайцев
        for (auto& rabbit : rabbits_) {
            rabbit.Move(n, m);
        }

        if (search_) {
            // Движение волков с поиском зайцев
            for (auto& wolf : wolves_) {
                wolf.MoveWithSearch(n, m, rabbits_);
            }
        }
        else {
            // Движение волков без поиска
            for (auto& wolf : wolves_) {
                wolf.MoveWithoutSearch(n, m);
            }
        }


        // Размножение зайцев
        std::vector<Rabbit> newRabbits;
        for (auto& rabbit : rabbits_) {
            if (rabbit.getAge() >= max_age_of_rabbits) {
                rabbit.resetAge();
                if (rabbits_.size() + newRabbits.size() <max_rabbits_count) {
                    newRabbits.emplace_back(rabbit.getX(), rabbit.getY());
                }
            }
        }
        rabbits_.insert(rabbits_.end(), newRabbits.begin(), newRabbits.end());

        // Взаимодействие волков и зайцев
        std::vector<Wolf> newWolves;
        std::vector<Wolf> survivingWolves;

        for (auto& wolf : wolves_) {
            bool ateRabbit = false;
            auto it = rabbits_.begin();
            while (it != rabbits_.end()) {
                if (wolf.getPosition() == it->getPosition()) {
                    ateRabbit = true;
                    wolf.resetAge();
                    if(wolves_.size() + newWolves.size() < max_wolf_count){
                        newWolves.emplace_back(wolf.getX(), wolf.getY());
                    }
                    it = rabbits_.erase(it);
                    break;
                }
                else {
                    ++it;
                }
            }

            if (ateRabbit || wolf.getAge() < max_age_of_wolfs) {
                survivingWolves.push_back(wolf);
            }
            // Волки с возрастом >= 9 и не съевшие зайца - погибают (не добавляются в survivingWolves)
        }

        wolves_ = survivingWolves;
        wolves_.insert(wolves_.end(), newWolves.begin(), newWolves.end());
    }

    void printStats() const {
        std::cout << "Зайцы: " << rabbits_.size() << ", Волки: " << wolves_.size() << std::endl;
    }

    bool isFinished() const {
        if (rabbits_.empty() && wolves_.empty()) {
            return true;
        }
        if (rabbits_.empty()) {
            std::cout << "Остались только Волки: " << wolves_.size() << std::endl;
            return true;
        }
        if (wolves_.empty()) {
            std::cout << "Остались только Зайцы: " << rabbits_.size() << std::endl;
            return true;
        }
        return false;
    }

    double calculateCorrelation() const {
        if (rabbits_.empty() || wolves_.empty()) return 0.0;

        // Рассчитываем средние количества животных на клетку
        double rabbitMean = static_cast<double>(rabbits_.size()) / (n * m);
        double wolfMean = static_cast<double>(wolves_.size()) / (n * m);

        // Создаем сетки для подсчета животных в каждой клетке
        std::vector<std::vector<int>> rabbitGrid(n, std::vector<int>(m, 0));
        std::vector<std::vector<int>> wolfGrid(n, std::vector<int>(m, 0));

        // Заполняем сетки
        for (const auto& rabbit : rabbits_) {
            int x = rabbit.getX() % n;
            int y = rabbit.getY() % m;
            if (x >= 0 && x < n && y >= 0 && y < m) {
                rabbitGrid[x][y]++;
            }
        }

        for (const auto& wolf : wolves_) {
            int x = wolf.getX() % n;
            int y = wolf.getY() % m;
            if (x >= 0 && x < n && y >= 0 && y < m) {
                wolfGrid[x][y]++;
            }
        }

        // Рассчитываем ковариацию
        double covariance = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                covariance += (rabbitGrid[i][j] - rabbitMean) * (wolfGrid[i][j] - wolfMean);
            }
        }
        covariance /= (n * m);

        // Рассчитываем стандартные отклонения
        double rabbitVariance = 0.0;
        double wolfVariance = 0.0;

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                rabbitVariance += std::pow(rabbitGrid[i][j] - rabbitMean, 2);
                wolfVariance += std::pow(wolfGrid[i][j] - wolfMean, 2);
            }
        }
        rabbitVariance /= (n * m);
        wolfVariance /= (n * m);

        double rabbitStd = std::sqrt(std::max(0.0, rabbitVariance));
        double wolfStd = std::sqrt(std::max(0.0, wolfVariance));

        if (rabbitStd == 0 || wolfStd == 0) return 0.0;

        return covariance / (rabbitStd * wolfStd);
    }

private:
    int n, m;
    std::vector<Rabbit> rabbits_;
    std::vector<Wolf> wolves_;
    bool search_ = false;
};

int main() {
    setlocale(LC_ALL, "Russian");

    const int N = size_of_iland;
    const int M = size_of_iland;


    //Поиск?
    Iland island(N, M, is_searching);



    // Начальная популяция
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, N - 1);

    for (int i = 0; i < start_count_of_rabbits; ++i) {
        island.addRabbit(dis(gen), dis(gen));
    }
    for (int i = 0; i < start_count_of_wolfs; ++i) {
        island.addWolf(dis(gen), dis(gen));
    }

    

    // Симуляция
    for (int step = 0; step < STEPS; ++step) {
        std::cout << "Шаг " << step + 1 << ": ";
        island.simulateStep();
        island.printStats();

        double correlation = island.calculateCorrelation();
        std::cout << "Корреляция: " << correlation << std::endl;

        if (island.isFinished()) {
            break;
        }

        std::cout << std::endl;
    }

    return 0;
}
