#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <list>
#include <memory>

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
const int GRID_CELL_SIZE = 10;

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

    int gridKey(int cellSize) const {
        return (x / cellSize) * 1000 + (y / cellSize);
    }
};

namespace std {
    template <>
    struct hash<Position> {
        size_t operator()(const Position& pos) const {
            return hash<int>()(pos.x) ^ (hash<int>()(pos.y) << 1);
        }
    };
}

class Animal {
protected:
    Position position_;
    int age;

public:
    Animal(int startX, int startY) : position_(startX, startY), age(0) {}
    Animal(Position startP) : position_(startP), age(0) {}
    virtual ~Animal() = default;

    virtual void Move(std::mt19937& gen, int width, int height) {
        std::uniform_int_distribution<> dis(-1, 1);
        position_.x += dis(gen);
        position_.y += dis(gen);

        position_.x = (position_.x + width) % width;
        position_.y = (position_.y + height) % height;

        age++;
    }

    int getX() const { return position_.x; }
    int getY() const { return position_.y; }
    Position getPosition() const { return position_; }
    int getAge() const { return age; }
    void resetAge() { age = 0; }
};

class Rabbit : public Animal {
public:
    Rabbit(int x, int y) : Animal(x, y) {}
    using Animal::Animal;
};

class Wolf : public Animal {
public:
    Wolf(int x, int y) : Animal(x, y) {}
    using Animal::Animal;

    void MoveWithSearch(std::mt19937& gen, int width, int height, const std::unordered_map<int, std::vector<Position>>& rabbitGrid, int cellSize) {
        int currentGridKey = position_.gridKey(cellSize);
        Position bestTarget;
        double bestDistance = std::numeric_limits<double>::max();
        bool found = false;

        int gridX = position_.x / cellSize;
        int gridY = position_.y / cellSize;

        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                int neighborKey = (gridX + dx) * 1000 + (gridY + dy);
                auto it = rabbitGrid.find(neighborKey);
                if (it != rabbitGrid.end()) {
                    for (const auto& rabbitPos : it->second) {
                        double dist = position_.distanceTo(rabbitPos);
                        if (dist < bestDistance) {
                            bestDistance = dist;
                            bestTarget = rabbitPos;
                            found = true;
                        }
                    }
                }
            }
        }

        if (found) {
            std::uniform_real_distribution<> prob(0.0, 1.0);

            if (prob(gen) < quality_of_search) {
                // Двигаемся к цели
                int dx = 0, dy = 0;
                if (bestTarget.x > position_.x) dx = 1;
                else if (bestTarget.x < position_.x) dx = -1;

                if (bestTarget.y > position_.y) dy = 1;
                else if (bestTarget.y < position_.y) dy = -1;

                position_.x += dx;
                position_.y += dy;
            }
            else {
                std::uniform_int_distribution<> dis(-1, 1);
                position_.x += dis(gen);
                position_.y += dis(gen);
            }
        }
        else {
            std::uniform_int_distribution<> dis(-1, 1);
            position_.x += dis(gen);
            position_.y += dis(gen);
        }

        position_.x = (position_.x + width) % width;
        position_.y = (position_.y + height) % height;
        age++;
    }
};

class Iland {
private:
    int n, m;
    std::vector<Rabbit> rabbits_;
    std::vector<Wolf> wolves_;
    bool search_;

    std::unordered_map<int, std::vector<Position>> rabbitGrid_;
    int cellSize_;

    std::mt19937 gen_;

    void updateRabbitGrid() {
        rabbitGrid_.clear();
        for (const auto& rabbit : rabbits_) {
            Position pos = rabbit.getPosition();
            rabbitGrid_[pos.gridKey(cellSize_)].push_back(pos);
        }
    }

public:
    Iland(int width, int height, bool search)
        : n(width), m(height), search_(search), cellSize_(GRID_CELL_SIZE), gen_(std::random_device{}()) {
    }

    void addRabbit(int x, int y) {
        x = (x + n) % n;
        y = (y + m) % m;
        rabbits_.emplace_back(x, y);
    }

    void addWolf(int x, int y) {
        wolves_.emplace_back(x, y);
    }

    void simulateStep() {
        for (auto& rabbit : rabbits_) {
            rabbit.Move(gen_, n, m);
        }
        updateRabbitGrid();

        if (search_) {
            for (auto& wolf : wolves_) {
                wolf.MoveWithSearch(gen_, n, m, rabbitGrid_, cellSize_);
            }
        }
        else {
            for (auto& wolf : wolves_) {
                wolf.Move(gen_, n, m);
            }
        }

        size_t currentRabbits = rabbits_.size();
        if (currentRabbits < max_rabbits_count) {
            std::vector<Rabbit> newRabbits;
            newRabbits.reserve(currentRabbits / 2);

            for (auto& rabbit : rabbits_) {
                if (rabbit.getAge() >= max_age_of_rabbits) {
                    rabbit.resetAge();
                    if (rabbits_.size() + newRabbits.size() < max_rabbits_count) {
                        Position pos = rabbit.getPosition();
                        newRabbits.emplace_back(pos.x, pos.y);
                    }
                }
            }

            if (!newRabbits.empty()) {
                rabbits_.insert(rabbits_.end(), newRabbits.begin(), newRabbits.end());
                updateRabbitGrid();
            }
        }

        std::vector<Wolf> newWolves;
        newWolves.reserve(wolves_.size() / 2);

        std::vector<Wolf> survivingWolves;
        survivingWolves.reserve(wolves_.size());

        std::unordered_map<Position, size_t, std::hash<Position>> rabbitPositions;
        for (size_t i = 0; i < rabbits_.size(); ++i) {
            rabbitPositions[rabbits_[i].getPosition()] = i;
        }

        for (auto& wolf : wolves_) {
            Position wolfPos = wolf.getPosition();
            auto it = rabbitPositions.find(wolfPos);

            if (it != rabbitPositions.end()) {
                // Волк съел зайца
                wolf.resetAge();
                if (wolves_.size() + newWolves.size() < max_wolf_count) {
                    newWolves.emplace_back(wolf.getX(), wolf.getY());
                }

                size_t rabbitIndex = it->second;
                if (rabbitIndex != rabbits_.size() - 1) {
                    std::swap(rabbits_[rabbitIndex], rabbits_.back());
                    rabbitPositions[rabbits_[rabbitIndex].getPosition()] = rabbitIndex;
                }
                rabbits_.pop_back();
                rabbitPositions.erase(it);

                updateRabbitGrid();
            }

            if (wolf.getAge() < max_age_of_wolfs) {
                survivingWolves.push_back(std::move(wolf));
            }
        }

        wolves_ = std::move(survivingWolves);
        if (!newWolves.empty()) {
            wolves_.insert(wolves_.end(), newWolves.begin(), newWolves.end());
        }
    }

    void printStats() const {
        std::cout << "Зайцы: " << rabbits_.size() << ", Волки: " << wolves_.size() << std::endl;
    }

    bool isFinished() const {
        if (rabbits_.empty() && wolves_.empty()) return true;
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

        // Используем уже рассчитанную сетку для зайцев
        std::unordered_map<int, int> wolfGrid;
        for (const auto& wolf : wolves_) {
            wolfGrid[wolf.getPosition().gridKey(cellSize_)]++;
        }

        double rabbitMean = static_cast<double>(rabbits_.size()) / (n * m);
        double wolfMean = static_cast<double>(wolves_.size()) / (n * m);

        double covariance = 0.0;
        double rabbitVariance = 0.0;
        double wolfVariance = 0.0;

        // Проходим по всем ячейкам сетки
        for (const auto& rabbitCell : rabbitGrid_) {
            int wolfCount = wolfGrid[rabbitCell.first];
            double rabbitCellCount = rabbitCell.second.size() / static_cast<double>(cellSize_ * cellSize_);
            double wolfCellCount = wolfCount / static_cast<double>(cellSize_ * cellSize_);

            covariance += (rabbitCellCount - rabbitMean) * (wolfCellCount - wolfMean);
            rabbitVariance += std::pow(rabbitCellCount - rabbitMean, 2);
            wolfVariance += std::pow(wolfCellCount - wolfMean, 2);
        }

        int totalCells = (n / cellSize_) * (m / cellSize_);
        if (totalCells > 0) {
            covariance /= totalCells;
            rabbitVariance /= totalCells;
            wolfVariance /= totalCells;
        }

        double rabbitStd = std::sqrt(std::max(0.0, rabbitVariance));
        double wolfStd = std::sqrt(std::max(0.0, wolfVariance));

        if (rabbitStd == 0 || wolfStd == 0) return 0.0;

        return covariance / (rabbitStd * wolfStd);
    }
};

int main() {
    setlocale(LC_ALL, "Russian");

    const int N = size_of_iland;
    const int M = size_of_iland;

    Iland island(N, M, is_searching);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, N - 1);

    for (int i = 0; i < start_count_of_rabbits; ++i) {
        island.addRabbit(dis(gen), dis(gen));
    }
    for (int i = 0; i < start_count_of_wolfs; ++i) {
        island.addWolf(dis(gen), dis(gen));
    }

    for (int step = 0; step < STEPS; ++step) {
        std::cout << "Шаг " << step + 1 << ": ";
        island.simulateStep();
        island.printStats();

        double correlation = island.calculateCorrelation();
        std::cout << "Корреляция: " << correlation << std::endl;

        if (island.isFinished()) {
            break;
        }
    }

    return 0;
}
